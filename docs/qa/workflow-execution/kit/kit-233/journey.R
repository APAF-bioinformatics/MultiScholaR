# KIT-233 (r1.2): the MultiScholaR application driven in a headless browser through the default DIA proteomics
# journey, on a GitHub-hosted runner (run-tests scripts job). Records only: nothing under R/ is changed, and what the
# journey finds is reported, never fixed.
#
# As a script, from the package root:
#   Rscript docs/qa/workflow-execution/kit/kit-233/journey.R lifecycle <commit under test>   run L (L-0 to L-5)
#   Rscript docs/qa/workflow-execution/kit/kit-233/journey.R journey <commit under test>     run J (PXD055981 s500)
#   Rscript docs/qa/workflow-execution/kit/kit-233/journey.R userlife <commit under test>    L-u (user mode, observation)
#   Rscript docs/qa/workflow-execution/kit/kit-233/journey.R journey-inner <phase file>      the journey phase (5.1)
# Sourced (`source(<this file>)`), it only defines the functions and the route table (KIT-237 need N1). The launch is a
# parameter (`k233LaunchSpec`: the package tree's entry.R, or an installed application folder, with the network method
# recorded; N2), and the steps come from a route table (`K233_ROUTE_DIA`; N3).
# Environment: K233_LAUNCH_FORM = "e" (default: `Rscript --vanilla -e ...`) or "wrapper" (5.1's recorded fallback).

K233_OUT <- "k233-out"
K233_RECORDS <- "docs/qa/workflow-execution/kit/kit-233"
K233_NS <- "proteomics_workflow-"
K233_PROBE <- '{"launcher":"process","reason":"KIT-233 hosted runner: no user D-Bus"}'
K233_TEST_OPTIONS <- "options(shiny.testmode = TRUE, shiny.maxRequestSize = 2 * 1024^3)"
K233_TREE_ENTRY <- "packaging/electron/entry.R"
K233_NO_BROWSER <- "/nonexistent/k233-no-browser-for-the-application"
K233_PXD <- "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/12/PXD055981/"
K233_INPUTS <- data.frame(
    role = c("diann_report", "fasta_human", "fasta_yeast", "fasta_ecoli"),
    file = c(
        "ultra-lfq-all-legacy.parquet",
        "human_canonical_uniprotkb_proteome_UP000005640_2023_12_16.fasta",
        "yeast_canonical_uniprotkb_proteome_UP000002311_2023_12_16.fasta",
        "Ecoli_canonical_uniprotkb_proteome_UP000000625_2023_12_16.fasta"
    ),
    bytes = c(1229341213, 39640337, 3853617, 1890883),
    stringsAsFactors = FALSE
)
K233_S500 <- ".*_500_DIA_15_(H1E-[BCDE])_.*"
K233_MIXES <- c("H1E-B", "H1E-C", "H1E-D", "H1E-E")
K233_PANDOC <- "3.11" # distributable.yml:72, the bundle's pinned pandoc (packaging/linux/build.sh:95-101)

# ---------------------------------------------------------------------------------------------------------------
# The DIA route (N3). Titles and tab ids: R/utils_omic_gui_config.R:2236-2258 (the application's own table; the
# driver reads the live table first and uses this one only if that read fails). Apply buttons and R6 states:
# tests/testthat/helper-e2e-browser.R:1263-1287, :1311-1325. Result and plot outputs: the sub-tab UIs
# (R/mod_prot_qc_peptide_*.R, R/mod_prot_qc_protein_*.R). Permissive values: helper-e2e-browser.R:1266-1282,
# :1315-1324, applied only by L-1's rule (a default refused on the 8-run fixture).
# ---------------------------------------------------------------------------------------------------------------
K233_ROUTE_DIA <- list(
    id = "DIA",
    steps = list(
        list(id = "S1", name = "home and project", fn = "k233StepHome"),
        list(id = "S2", name = "import", fn = "k233StepImport"),
        list(id = "S3", name = "design", fn = "k233StepDesign"),
        list(id = "S4", name = "peptide QC", fn = "k233StepQc", group = "peptide_qc", title = "Peptide QC",
            tabset = "peptide_qc-peptide_filter_tabs"),
        list(id = "S5", name = "protein QC", fn = "k233StepQc", group = "protein_qc", title = "Protein QC",
            tabset = "protein_qc-protein_filter_tabs"),
        list(id = "S6", name = "normalisation", fn = "k233StepNormalisation"),
        list(id = "S7", name = "export", fn = "k233StepExport"),
        list(id = "S8", name = "differential abundance", fn = "k233StepDa"),
        list(id = "S9", name = "enrichment", fn = "k233StepEnrichment", continue_on_refusal = TRUE),
        list(id = "S10", name = "summary", fn = "k233StepSummary")
    ),
    qc_tabset = "qc_tabs_lfq",
    qc = list(
        peptide_qc = list(
            list(id = "qvalue_filter", title = "Q-Value Filter", button = "apply_qvalue_filter",
                state = "qvalue_filtered", outputs = c("qvalue_results", "qvalue_plot")),
            list(id = "protein_peptide_filter", title = "Protein Evidence", button = "apply_protein_peptide_filter",
                state = "protein_peptide_filtered", outputs = c("protein_peptida_results", "protein_peptide_plot"),
                permissive = list(min_peptides_per_protein = 1, min_peptidoforms_per_protein = 1)),
            list(id = "rollup", title = "Precursor Rollup", button = "apply_rollup", state = "precursor_rollup",
                outputs = c("rollup_results", "rollup_plot")),
            list(id = "intensity_filter", title = "Intensity Filter", button = "apply_intensity_filter",
                state = "intensity_filtered", outputs = c("intensity_results", "intensity_plot"),
                permissive = list(use_strict_mode = FALSE, min_reps_per_group = 1, min_groups = 1,
                    intensity_cutoff_percentile = 0.1)),
            list(id = "sample_filter", title = "Sample Quality", button = "apply_sample_filter",
                state = "sample_filtered", outputs = c("sample_results", "sample_plot"),
                permissive = list(min_peptides_per_sample = 1)),
            list(id = "replicate_filter", title = "Replicate Filter", button = "apply_replicate_filter",
                state = "replicate_filtered", outputs = c("replicate_results", "replicate_plot")),
            list(id = "imputation", title = "Imputation", button = "apply_imputation", state = "imputed",
                outputs = c("imputation_results", "imputation_plot"),
                permissive = list(proportion_missing_values = 1.0))
        ),
        protein_qc = list(
            list(id = "rollup", title = "IQ Protein Rollup", button = "apply_iq_rollup", state = "protein_s4_created",
                outputs = c("iq_rollup_results", "iq_rollup_plot")),
            list(id = "cleanup", title = "Accession Cleanup", button = "apply_accession_cleanup",
                state = "protein_accession_cleaned", outputs = c("accession_cleanup_results", "accession_cleanup_plot")),
            list(id = "intensity_filter", title = "Protein Intensity Filter",
                button = "apply_protein_intensity_filter", state = "protein_intensity_filtered",
                outputs = c("protein_intensity_filter_results", "protein_intensity_filter_plot"),
                permissive = list(use_strict_mode = FALSE, min_reps_per_group = 1, min_groups = 1,
                    proteins_intensity_cutoff_percentile = 0.1)),
            list(id = "duplicate_removal", title = "Duplicate Removal", button = "apply_duplicate_removal",
                state = "duplicates_removed", outputs = c("duplicate_removal_results", "duplicate_removal_plot")),
            list(id = "replicate_filter", title = "Protein Replicate Filter",
                button = "apply_protein_replicate_filter", state = "protein_replicate_filtered",
                outputs = c("protein_replicate_filter_results", "protein_replicate_filter_plot"),
                permissive = list(parallel_cores = 1))
        )
    )
)

# ---------------------------------------------------------------------------------------------------------------
# Small utilities.
# ---------------------------------------------------------------------------------------------------------------
k233Clean <- function(x) {
    x <- if (is.null(x) || length(x) == 0L) "" else paste(as.character(x), collapse = " | ")
    gsub("[\t\r\n]+", " | ", x)
}

k233Bind <- function(rows) {
    cols <- unique(unlist(lapply(rows, names)))
    columns <- lapply(stats::setNames(cols, cols), \(k) vapply(rows, \(r) k233Clean(r[[k]]), character(1)))
    data.frame(columns, check.names = FALSE, stringsAsFactors = FALSE)
}

k233WriteTsv <- function(rows, path) {
    if (is.data.frame(rows)) {
        rows[] <- lapply(rows, \(col) gsub("[\t\r\n]+", " | ", as.character(col)))
        utils::write.table(rows, path, sep = "\t", quote = FALSE, row.names = FALSE)
    } else if (length(rows) == 0L) {
        writeLines("(no rows)", path)
    } else {
        utils::write.table(k233Bind(rows), path, sep = "\t", quote = FALSE, row.names = FALSE)
    }
    invisible(path)
}

k233Sha <- function(path) {
    if (length(path) == 1L && file.exists(path)) digest::digest(file = path, algo = "sha256") else NA_character_
}

k233Utc <- function(time = Sys.time()) format(time, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

k233Seconds <- function(since) as.numeric(difftime(Sys.time(), since, units = "secs"))

k233ReadLog <- function(path) if (file.exists(path)) readLines(path, warn = FALSE) else character()

k233Norm <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))

k233PngSize <- function(path) {
    if (!file.exists(path)) return(c(NA_real_, NA_real_))
    b <- readBin(path, "raw", 24L)
    if (length(b) < 24L) return(c(NA_real_, NA_real_))
    big <- function(x) sum(as.numeric(as.integer(x)) * 256^(3:0))
    c(big(b[17:20]), big(b[21:24]))
}

k233Bounds <- function(run_mode) {
    if (identical(run_mode, "journey")) {
        return(list(deadline_min = 155, hard_s = 35 * 60, finding_s = 600, poll = 1))
    }
    list(deadline_min = 40, hard_s = 5 * 60, finding_s = 600, poll = 0.5)
}

#' The launch (N2): `form = "tree"` starts the package tree's `packaging/electron/entry.R` (r1.2 5.1); `form =
#' "installed"` starts `entry.R` under an installed application folder, as `val-parity.R:64-65` finds it. `network`
#' is the method recorded and, for "namespace", the one the outer phase applies (5.1); "proxy" and "firewall" are
#' recorded as given (a caller that closes the network itself names its method).
k233LaunchSpec <- function(form = "tree", app_dir = NULL, network = "namespace") {
    list(form = form, app_dir = app_dir, network = network, e_form = Sys.getenv("K233_LAUNCH_FORM", "e"))
}

k233EntryPath <- function(launch) {
    if (!identical(launch$form, "installed")) return(K233_TREE_ENTRY)
    entry <- file.path(launch$app_dir, c("resources", "app/resources"), "entry.R")
    entry <- entry[file.exists(entry)]
    if (!length(entry)) stop("no entry.R under the installed application folder ", launch$app_dir)
    normalizePath(entry[[1L]], winslash = "/")
}

# ---------------------------------------------------------------------------------------------------------------
# Outer phase: environment, the commit check, J's data phase, then the journey phase in a network namespace.
# ---------------------------------------------------------------------------------------------------------------
k233Main <- function(args) {
    if (length(args) < 2L) {
        stop("usage: journey.R <lifecycle|journey|userlife> <commit under test> | journey-inner <phase file>")
    }
    if (identical(args[[1L]], "journey-inner")) return(invisible(k233Inner(args[[2L]])))
    if (!args[[1L]] %in% c("lifecycle", "journey", "userlife")) stop("unknown run mode ", args[[1L]])
    k233Outer(args[[1L]], args[[2L]], k233LaunchSpec())
}

k233ScriptPath <- function() {
    file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    if (length(file)) normalizePath(file[[1L]]) else normalizePath(file.path(K233_RECORDS, "journey.R"))
}

k233OuterEnvironment <- function(run_mode, commit, pkg, launch) {
    git <- function(...) suppressWarnings(system2("git", c("-C", shQuote(pkg), ...), stdout = TRUE, stderr = TRUE))
    head <- git("rev-parse", "HEAD")
    changed <- git("diff", "--name-only", commit, "HEAD")
    diff_ok <- is.null(attr(changed, "status")) && all(startsWith(changed, paste0(K233_RECORDS, "/")))
    chrome <- tryCatch(chromote::find_chrome(), error = \(e) "")
    chrome_version <- if (nzchar(chrome)) {
        tryCatch(suppressWarnings(system2(chrome, "--version", stdout = TRUE, stderr = TRUE)), error = \(e) "?")
    } else {
        "not found"
    }
    version <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = \(e) "not installed")
    cpu <- grep("^model name", tryCatch(readLines("/proc/cpuinfo"), error = \(e) character()), value = TRUE)
    mem <- grep("^MemTotal", tryCatch(readLines("/proc/meminfo"), error = \(e) character()), value = TRUE)
    rows <- list(
        list(key = "run_mode", value = run_mode),
        list(key = "mode", value = if (identical(run_mode, "userlife")) "user" else "test"),
        list(key = "commit_under_test", value = commit),
        list(key = "records_commit", value = head),
        list(key = "diff_commit_under_test_vs_records", value = paste(changed, collapse = " ")),
        list(key = "diff_only_records_folder", value = diff_ok),
        list(key = "launch_form", value = launch$form),
        list(key = "launch_e_form", value = launch$e_form),
        list(key = "network_requested", value = launch$network),
        list(key = "r_version", value = R.version.string),
        list(key = "effective_user_outer", value = Sys.info()[["effective_user"]]),
        list(key = "blas", value = extSoftVersion()[["BLAS"]]),
        list(key = "lapack", value = La_library()),
        list(key = "threads_outer", value = paste0("OMP_NUM_THREADS=", Sys.getenv("OMP_NUM_THREADS"),
            " OPENBLAS_NUM_THREADS=", Sys.getenv("OPENBLAS_NUM_THREADS"))),
        list(key = "threads_application", value = "OMP/OPENBLAS/MKL_NUM_THREADS=1 (5.1); workers set 1 themselves"),
        list(key = "chrome", value = chrome),
        list(key = "chrome_version", value = chrome_version),
        list(key = "shinytest2", value = version("shinytest2")),
        list(key = "chromote", value = version("chromote")),
        list(key = "processx", value = version("processx")),
        list(key = "MultiScholaR", value = version("MultiScholaR")),
        list(key = "cpu", value = if (length(cpu)) sub("^model name\\s*:\\s*", "", cpu[[1L]]) else "?"),
        list(key = "cores", value = parallel::detectCores()),
        list(key = "memory", value = if (length(mem)) mem[[1L]] else "?"),
        list(key = "github_run_id", value = Sys.getenv("GITHUB_RUN_ID")),
        list(key = "started_utc", value = k233Utc())
    )
    list(rows = rows, commit_ok = diff_ok, head = head)
}

k233Outer <- function(run_mode, commit, launch) {
    started <- Sys.time()
    pkg <- normalizePath(".")
    out <- file.path(pkg, K233_OUT)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    work <- file.path(tempdir(), "k233-work")
    dir.create(work, recursive = TRUE, showWarnings = FALSE)
    env <- k233OuterEnvironment(run_mode, commit, pkg, launch)
    k233WriteTsv(env$rows, file.path(out, "environment.tsv"))
    if (!isTRUE(env$commit_ok)) {
        writeLines("stopped: the records commit differs from the commit under test outside the records folder (8)",
            file.path(out, "stop.txt"))
        return(k233Finish(out))
    }
    inputs <- list(id = "prot_dia fixture", sha = k233Sha(file.path(pkg, "tests/testdata/e2e/prot_dia/report.tsv")))
    if (identical(run_mode, "journey")) {
        data <- tryCatch(k233DataPhase(pkg, out, work), error = \(e) list(ok = FALSE,
            stop = paste("data phase error:", conditionMessage(e))))
        if (!isTRUE(data$ok)) {
            writeLines(data$stop, file.path(out, "stop.txt"))
            return(k233Finish(out))
        }
        inputs <- data$inputs
    }
    if (identical(run_mode, "userlife")) inputs$home_dir <- k233PlaceUserInputs(pkg)
    pandoc <- k233Pandoc(work)
    phase <- list(
        run_mode = run_mode, started = as.numeric(started), bounds = k233Bounds(run_mode), pkg = pkg, out = out,
        work = work, script = k233ScriptPath(), commit = commit, records_commit = env$head, inputs = inputs,
        pandoc = pandoc, launch = launch, user_outer = Sys.info()[["effective_user"]], libpaths = .libPaths(),
        locale = Sys.getlocale("LC_CTYPE"),
        env = c(
            R_LIBS = Sys.getenv("R_LIBS"), R_LIBS_USER = Sys.getenv("R_LIBS_USER"),
            R_LIBS_SITE = Sys.getenv("R_LIBS_SITE"), LANG = Sys.getenv("LANG", "C.UTF-8"),
            GITHUB_RUN_ID = Sys.getenv("GITHUB_RUN_ID"), TZ = "UTC",
            OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
            # shinytest2::AppDriver skips unless NOT_CRAN is true (run L 45523435786: "Reason: On CRAN").
            NOT_CRAN = "true"
        )
    )
    phase$env <- phase$env[nzchar(phase$env)]
    phase_file <- file.path(work, "k233-phase.rds")
    saveRDS(phase, phase_file)
    k233RunInner(phase, phase_file)
    inner_env <- file.path(out, "environment-inner.tsv")
    if (file.exists(inner_env)) {
        inner <- utils::read.delim(inner_env, stringsAsFactors = FALSE)
        user_inner <- inner$value[inner$key == "effective_user_inner"]
        rows <- c(env$rows, list(
            list(key = "effective_user_inner", value = user_inner),
            list(key = "effective_user_equal", value = identical(user_inner, phase$user_outer)),
            list(key = "pandoc", value = paste(pandoc$source, pandoc$dir))
        ))
        k233WriteTsv(rows, file.path(out, "environment.tsv"))
    }
    k233Finish(out)
}

k233RunInner <- function(phase, phase_file) {
    marker <- file.path(phase$out, "inner-started.txt")
    console <- file.path(phase$out, "inner-console.log")
    network <- phase$launch$network
    if (identical(network, "namespace")) {
        command <- paste(
            'ip link set lo up && exec sudo -u "$0" env PATH="$1" HOME="$2"',
            "http_proxy=http://127.0.0.1:9 https_proxy=http://127.0.0.1:9 no_proxy=127.0.0.1,localhost",
            'Rscript "$3" journey-inner "$4"'
        )
        user <- system2("id", "-un", stdout = TRUE)
        status <- system2("sudo", c("unshare", "--net", "--", "bash", "-c", shQuote(command), shQuote(user),
            shQuote(Sys.getenv("PATH")), shQuote(Sys.getenv("HOME")), shQuote(phase$script), shQuote(phase_file)),
            stdout = console, stderr = console)
        writeLines(sprintf("namespace command exit status %s", status), file.path(phase$out, "namespace-status.txt"))
        if (file.exists(marker)) return(invisible(status))
        network <- "proxy"
    }
    if (identical(network, "proxy")) {
        phase$launch$network <- "proxy"
        saveRDS(phase, phase_file)
        rscript <- file.path(R.home("bin"), "Rscript")
        status <- system2(rscript, c(shQuote(phase$script), "journey-inner", shQuote(phase_file)),
            stdout = file.path(phase$out, "inner-console-proxy.log"), stderr = file.path(phase$out, "inner-console-proxy.log"),
            env = c("http_proxy=http://127.0.0.1:9", "https_proxy=http://127.0.0.1:9", "no_proxy=127.0.0.1,localhost"))
        return(invisible(status))
    }
    rscript <- file.path(R.home("bin"), "Rscript")
    invisible(system2(rscript, c(shQuote(phase$script), "journey-inner", shQuote(phase_file)),
        stdout = console, stderr = console))
}

k233Finish <- function(out) {
    files <- setdiff(list.files(out, recursive = TRUE, all.files = TRUE), "index.tsv")
    index <- data.frame(path = files, bytes = file.size(file.path(out, files)),
        sha256 = vapply(file.path(out, files), k233Sha, character(1)), stringsAsFactors = FALSE)
    k233WriteTsv(index, file.path(out, "index.tsv"))
    for (f in c("stop.txt", "verdict.txt", "inner-fatal.txt")) {
        if (file.exists(file.path(out, f))) cat(f, ":", readLines(file.path(out, f)), sep = "\n")
    }
    if (file.exists(file.path(out, "steps.tsv"))) {
        steps <- utils::read.delim(file.path(out, "steps.tsv"), stringsAsFactors = FALSE)
        if (all(c("session", "step", "outcome", "seconds") %in% names(steps))) {
            print(steps[, c("session", "step", "outcome", "seconds")], row.names = FALSE)
        }
    }
    invisible(out)
}

k233Download <- function(url, dest) {
    started <- Sys.time()
    status <- system2("curl", c("-fsSL", "--retry", "5", "--retry-delay", "15", "-o", shQuote(dest), shQuote(url)))
    list(status = status, seconds = round(k233Seconds(started)), bytes = if (file.exists(dest)) file.size(dest) else NA)
}

#' J's data phase (7.1): downloads, sizes against PRIDE's listing, SHA-256, the s500 subset (stop rules 2 and 2b) and
#' the concatenated FASTA. Only the 500 SPD run names are written; no other series is subset or shown (stop rule 6).
k233DataPhase <- function(pkg, out, work) {
    dir <- file.path(work, "pxd055981")
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    rows <- list()
    for (i in seq_len(nrow(K233_INPUTS))) {
        dest <- file.path(dir, K233_INPUTS$file[[i]])
        url <- paste0(K233_PXD, K233_INPUTS$file[[i]])
        got <- k233Download(url, dest)
        rows[[i]] <- list(role = K233_INPUTS$role[[i]], file = K233_INPUTS$file[[i]], url = url,
            bytes_pride = K233_INPUTS$bytes[[i]], bytes = got$bytes, size_ok = isTRUE(got$bytes == K233_INPUTS$bytes[[i]]),
            sha256 = k233Sha(dest), seconds = got$seconds, rows = "", note = sprintf("curl exit %s", got$status))
    }
    k233WriteTsv(rows, file.path(out, "inputs.tsv"))
    if (!all(vapply(rows, \(r) isTRUE(r$size_ok), logical(1)))) {
        return(list(ok = FALSE, stop = "stop rule 2: a download's size differs from PRIDE's listing (inputs.tsv)"))
    }
    frozen <- file.path(pkg, "docs/qa/workflow-execution/kit/kit-160/inputs.tsv")
    r0 <- "no KIT-160 R0 record at the commit under test; hashes recorded here for R0 to compare"
    if (file.exists(frozen)) {
        record <- utils::read.delim(frozen, stringsAsFactors = FALSE)
        if (all(c("file", "sha256") %in% names(record))) {
            for (r in rows) {
                want <- record$sha256[record$file == r$file]
                if (length(want) && nzchar(want[[1L]]) && !identical(want[[1L]], r$sha256)) {
                    return(list(ok = FALSE, stop = paste("stop rule 2: SHA-256 differs from KIT-160's frozen", r$file)))
                }
            }
            r0 <- "checked against KIT-160's frozen inputs.tsv: equal where frozen"
        }
    }
    report <- file.path(dir, K233_INPUTS$file[[1L]])
    ds <- arrow::open_dataset(report)
    counts <- tryCatch(ds |> dplyr::count(Run) |> dplyr::collect() |> as.data.frame(),
        error = \(e) as.data.frame(table(Run = arrow::read_parquet(report, col_select = "Run")$Run), responseName = "n"))
    runs <- as.character(counts$Run)
    s500 <- sort(grep(K233_S500, runs, value = TRUE))
    mix <- sub(K233_S500, "\\1", s500)
    found <- table(factor(mix, levels = K233_MIXES))
    if (length(s500) != 16L || any(found != 4L)) {
        names500 <- sort(grep("_500_", runs, value = TRUE))
        k233WriteTsv(data.frame(run_500_spd = names500), file.path(out, "runs-500spd-found.tsv"))
        return(list(ok = FALSE, stop = sprintf(
            "stop rule 2b: the s500 pattern gave %d runs (%s), not 16 with 4 per mix; 500 SPD names in runs-500spd-found.tsv",
            length(s500), paste(names(found), found, sep = "=", collapse = ", "))))
    }
    subset_path <- file.path(dir, "pxd055981-s500-15ng.parquet")
    subset <- ds |> dplyr::filter(Run %in% s500) |> dplyr::collect()
    arrow::write_parquet(subset, subset_path)
    n_rows <- nrow(subset)
    rm(subset)
    gc()
    fasta <- file.path(dir, "pxd055981-human-yeast-ecoli.fasta")
    k233Concatenate(file.path(dir, K233_INPUTS$file[2:4]), fasta)
    run_rows <- counts[counts$Run %in% s500, , drop = FALSE]
    k233WriteTsv(data.frame(run = run_rows$Run, mix = sub(K233_S500, "\\1", run_rows$Run), rows = run_rows$n),
        file.path(out, "inputs-runs.tsv"))
    rows <- c(rows, list(
        list(role = "subset_s500", file = basename(subset_path), url = "derived on the runner (7.1)", bytes_pride = "",
            bytes = file.size(subset_path), size_ok = "", sha256 = k233Sha(subset_path), seconds = "", rows = n_rows,
            note = paste("16 runs, 4 per mix; run list in inputs-runs.tsv;", r0)),
        list(role = "fasta_concatenated", file = basename(fasta), url = "human, yeast, E. coli in that order (KIT-160 4.1)",
            bytes_pride = "", bytes = file.size(fasta), size_ok = "", sha256 = k233Sha(fasta), seconds = "", rows = "",
            note = "a newline is added between files that do not end with one")
    ))
    k233WriteTsv(rows, file.path(out, "inputs.tsv"))
    list(ok = TRUE, inputs = list(id = "PXD055981 s500", sha = k233Sha(subset_path), report = subset_path,
        fasta = fasta, runs = s500, mix = mix))
}

k233Concatenate <- function(files, dest) {
    con <- file(dest, "wb")
    on.exit(close(con))
    for (f in files) {
        bytes <- readBin(f, "raw", file.size(f))
        writeBin(bytes, con)
        if (length(bytes) && bytes[[length(bytes)]] != as.raw(10L)) writeBin(as.raw(10L), con)
    }
    invisible(dest)
}

#' L-u's inputs on the runner's disk, where the real file chooser reaches them ($HOME/k233-input).
k233PlaceUserInputs <- function(pkg) {
    home <- file.path(Sys.getenv("HOME"), "k233-input")
    dir.create(home, recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(pkg, "tests/testdata/e2e/prot_dia", c("report.tsv", "proteins.fasta")), home, overwrite = TRUE)
    home
}

#' Pandoc for the report as the product has it: the runner's when rmarkdown finds one, otherwise the bundle's pinned
#' version given to the application as RSTUDIO_PANDOC, as the shell does (packaging/electron/main.js:150).
k233Pandoc <- function(work) {
    if (isTRUE(tryCatch(rmarkdown::pandoc_available(), error = \(e) FALSE))) {
        return(list(dir = "", source = paste("runner pandoc", tryCatch(as.character(rmarkdown::pandoc_version()),
            error = \(e) "?"))))
    }
    url <- sprintf("https://github.com/jgm/pandoc/releases/download/%s/pandoc-%s-linux-amd64.tar.gz", K233_PANDOC,
        K233_PANDOC)
    tar <- file.path(work, "pandoc.tar.gz")
    got <- k233Download(url, tar)
    if (!identical(got$status, 0L)) return(list(dir = "", source = "none: the runner has no pandoc and the download failed"))
    utils::untar(tar, exdir = file.path(work, "pandoc"))
    dir <- file.path(work, "pandoc", paste0("pandoc-", K233_PANDOC), "bin")
    if (!file.exists(file.path(dir, "pandoc"))) return(list(dir = "", source = "none: pandoc archive without bin/pandoc"))
    list(dir = dir, source = sprintf("bundle-pinned pandoc %s via RSTUDIO_PANDOC (archive sha256 %s)", K233_PANDOC,
        k233Sha(tar)))
}

# ---------------------------------------------------------------------------------------------------------------
# Inner phase (inside the namespace): the sessions, the steps, the records.
# ---------------------------------------------------------------------------------------------------------------
k233Inner <- function(phase_file) {
    phase <- readRDS(phase_file)
    .libPaths(phase$libpaths)
    if (length(phase$env)) do.call(Sys.setenv, as.list(phase$env))
    Sys.unsetenv(c("MULTISCHOLAR_TEST_FEATURES", "MULTISCHOLAR_TEST_MODE", "CHROMOTE_CHROME"))
    invisible(tryCatch(Sys.setlocale("LC_ALL", phase$locale), warning = \(w) NULL, error = \(e) NULL))
    options(warn = 1, chromote.timeout = 60)
    writeLines(k233Utc(), file.path(phase$out, "inner-started.txt"))
    G <- k233Context(phase)
    on.exit({
        for (S in G$sessions) k233StopSession(S)
        k233WriteAll(G)
    }, add = TRUE)
    plan <- switch(phase$run_mode, lifecycle = k233Lifecycle, journey = k233Journey, userlife = k233UserLife)
    fatal <- function(e) {
        G$fatal <- paste(conditionMessage(e), "| class:", paste(class(e), collapse = ","), "| call:",
            paste(deparse(conditionCall(e)), collapse = " "))
        writeLines(G$fatal, file.path(G$out, "inner-fatal.txt"))
    }
    tryCatch(plan(G), error = fatal, skip = fatal)
    invisible(G)
}

k233Context <- function(phase) {
    G <- new.env(parent = emptyenv())
    G$phase <- phase
    G$out <- phase$out
    G$work <- phase$work
    G$pkg <- phase$pkg
    G$bounds <- phase$bounds
    G$deadline_at <- phase$started + phase$bounds$deadline_min * 60
    G$route <- K233_ROUTE_DIA
    for (k in c("steps", "findings", "numbers", "criteria", "sessions", "session_rows", "actions", "observations",
        "coverage")) {
        assign(k, list(), envir = G)
    }
    G$helpers <- new.env()
    sys.source(file.path(phase$pkg, "tests/testthat/helper-e2e-browser.R"), envir = G$helpers)
    links <- tryCatch(suppressWarnings(system2("ip", c("-o", "link"), stdout = TRUE, stderr = TRUE)),
        error = \(e) "ip not available")
    G$env <- list(
        list(key = "effective_user_inner", value = Sys.info()[["effective_user"]]),
        list(key = "network", value = phase$launch$network),
        list(key = "network_links", value = paste(sub("^[0-9]+: ([^:]+):.*$", "\\1", links), collapse = " ")),
        list(key = "network_closed_by", value = if (identical(phase$launch$network, "namespace")) {
            "network namespace with loopback only; proxies http://127.0.0.1:9"
        } else {
            paste(phase$launch$network, "(proxies http://127.0.0.1:9 where set)")
        }),
        list(key = "http_proxy", value = Sys.getenv("http_proxy")),
        list(key = "home", value = Sys.getenv("HOME")),
        list(key = "libpaths", value = paste(.libPaths(), collapse = ":")),
        list(key = "application_browser", value = paste("CHROMOTE_CHROME =", K233_NO_BROWSER, "(the application only)")),
        list(key = "gprofiler", value = if (identical(phase$run_mode, "userlife")) {
            "user mode: attempted and blocked (network closed)"
        } else {
            "test mode: not called (test-mode substitute)"
        })
    )
    G
}

k233Remaining <- function(G) G$deadline_at - as.numeric(Sys.time())

k233WriteAll <- function(G) {
    out <- G$out
    k233WriteTsv(G$steps, file.path(out, "steps.tsv"))
    k233WriteTsv(G$findings, file.path(out, "findings.tsv"))
    k233WriteTsv(G$numbers, file.path(out, "numbers.tsv"))
    k233WriteTsv(G$session_rows, file.path(out, "sessions.tsv"))
    k233WriteTsv(G$actions, file.path(out, "actions.tsv"))
    k233WriteTsv(G$observations, file.path(out, "observations.tsv"))
    k233WriteTsv(G$coverage, file.path(out, "coverage.tsv"))
    k233WriteTsv(G$env, file.path(out, "environment-inner.tsv"))
    if (length(G$criteria)) k233WriteTsv(G$criteria, file.path(out, "lifecycle.tsv"))
    k233ReviewIndex(G)
    invisible(G)
}

#' The reviewer's hand-off (7.3): each screenshot with its step name and time only.
k233ReviewIndex <- function(G) {
    rows <- list()
    for (r in G$steps) {
        shots <- grep("\\.png$", strsplit(k233Clean(r$files), " \\| ")[[1L]], value = TRUE)
        for (s in shots) {
            label <- if (startsWith(r$step, "S9")) {
                if (isTRUE(G$sessions[[r$session]]$test_mode)) "test-mode substitute" else "refused offline"
            } else {
                ""
            }
            rows[[length(rows) + 1L]] <- list(screenshot = sub(" .*$", "", s), step = r$step, name = r$name,
                seconds = r$seconds, label = label)
        }
    }
    k233WriteTsv(rows, file.path(G$out, "review-index.tsv"))
}

k233Criterion <- function(G, criterion, item, result, evidence = "") {
    if (length(result) != 1L) result <- FALSE
    value <- if (is.na(result)) "recorded" else if (isTRUE(result)) "pass" else "fail"
    G$criteria[[length(G$criteria) + 1L]] <- list(criterion = criterion, item = item, result = value,
        evidence = evidence)
    invisible(value)
}

k233Observation <- function(G, session, what, text) {
    G$observations[[length(G$observations) + 1L]] <- list(session = session, what = what, text = text,
        utc = k233Utc())
}

# --- sessions --------------------------------------------------------------------------------------------------

k233Launch <- function(G, label, test_mode = TRUE) {
    S <- new.env(parent = emptyenv())
    S$G <- G
    S$label <- label
    S$test_mode <- test_mode
    S$nn <- 0L
    S$hard_s <- G$bounds$hard_s
    S$poll <- G$bounds$poll
    S$hooks <- list()
    S$console_seen <- 0L
    S$port <- 18233L + length(G$sessions)
    S$run_folder <- file.path(G$work, paste0("run-", label))
    S$project_dir <- file.path(G$work, paste0("projects-", label))
    dir.create(S$run_folder, recursive = TRUE, showWarnings = FALSE)
    dir.create(S$project_dir, recursive = TRUE, showWarnings = FALSE)
    S$log <- file.path(G$out, sprintf("app-%s.log", label))
    G$sessions[[label]] <- S
    launch <- G$phase$launch
    entry <- k233EntryPath(launch)
    tail_args <- c(as.character(S$port), S$run_folder, as.character(Sys.getpid()))
    args <- if (!test_mode) {
        c("--vanilla", entry, tail_args)
    } else if (identical(launch$e_form, "wrapper")) {
        wrapper <- file.path(G$out, "k233-entry-wrapper.R")
        writeLines(c(K233_TEST_OPTIONS, sprintf('source("%s")', entry)), wrapper)
        c("--vanilla", wrapper, tail_args)
    } else {
        c("--vanilla", "-e", sprintf('%s; source("%s")', K233_TEST_OPTIONS, entry), tail_args)
    }
    env <- c("current", MSR_LAUNCHER_PROBE = K233_PROBE, OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1", CHROMOTE_CHROME = K233_NO_BROWSER)
    if (test_mode) env <- c(env, MULTISCHOLAR_TEST_MODE = "true")
    if (nzchar(G$phase$pandoc$dir)) env <- c(env, RSTUDIO_PANDOC = G$phase$pandoc$dir)
    S$command <- paste(c("Rscript", args), collapse = " ")
    started <- Sys.time()
    S$proc <- processx::process$new(file.path(R.home("bin"), "Rscript"), args, env = env, wd = G$pkg,
        stdout = S$log, stderr = "2>&1", cleanup = TRUE, cleanup_tree = TRUE)
    S$url <- sprintf("http://127.0.0.1:%d", S$port)
    S$answered <- FALSE
    for (i in seq_len(180L)) {
        if (!S$proc$is_alive()) break
        S$answered <- isTRUE(tryCatch(curl::curl_fetch_memory(S$url, handle = curl::new_handle(timeout = 5))$status_code ==
            200L, error = \(e) FALSE))
        if (S$answered) break
        Sys.sleep(1)
    }
    S$answer_s <- round(k233Seconds(started), 1)
    log <- k233ReadLog(S$log)
    S$launcher <- paste(grep("^MSR-LAUNCHER ", log, value = TRUE), collapse = " ")
    S$usage_error <- any(grepl("usage: Rscript --vanilla entry.R", log, fixed = TRUE))
    S$driver <- NULL
    S$attach_error <- ""
    if (S$answered) {
        refused <- function(e) {
            S$attach_error <- conditionMessage(e)
            NULL
        }
        S$driver <- tryCatch(
            shinytest2::AppDriver$new(S$url, load_timeout = 120000, timeout = 60000, width = 1440, height = 1100),
            error = refused, skip = refused
        )
    }
    S$recorder <- !is.null(S$driver) && isTRUE(k233Js(S, K233_RECORDER_JS))
    G$session_rows[[label]] <- list(session = label, mode = if (test_mode) "test" else "user", port = S$port,
        pid = S$proc$get_pid(), command = S$command, answered = S$answered, answer_seconds = S$answer_s,
        launcher = S$launcher, usage_error = S$usage_error, attached = !is.null(S$driver),
        attach_error = S$attach_error, recorder = S$recorder, started_utc = k233Utc(started))
    k233WriteAll(G)
    S
}

k233CopyLogs <- function(S) {
    files <- list.files(S$run_folder, pattern = "\\.log$", recursive = TRUE)
    for (f in files) {
        file.copy(file.path(S$run_folder, f), file.path(S$G$out, paste0(S$label, "-", gsub("/", "_", f))),
            overwrite = TRUE)
    }
    invisible(files)
}

k233StopSession <- function(S) {
    if (isTRUE(S$stopped)) return(invisible(NULL))
    S$stopped <- TRUE
    k233CopyLogs(S)
    if (!is.null(S$driver)) tryCatch(S$driver$stop(), error = \(e) NULL)
    tryCatch(S$proc$kill_tree(), error = \(e) tryCatch(S$proc$kill(), error = \(e) NULL))
    k233CopyLogs(S)
    G <- S$G
    if (!is.null(G$session_rows[[S$label]])) {
        G$session_rows[[S$label]]$exit_status <- tryCatch(S$proc$get_exit_status(), error = \(e) NA)
        G$session_rows[[S$label]]$stopped_utc <- k233Utc()
    }
    invisible(NULL)
}

# --- the page --------------------------------------------------------------------------------------------------

K233_RECORDER_JS <- r"---((function () {
  if (window.__k233) { return true; }
  var K = window.__k233 = { notes: [], modals: [] };
  function txt(e) { return ((e && (e.innerText || e.textContent)) || '').replace(/\s+/g, ' ').trim(); }
  function kind(c) {
    if (/shiny-notification-error/.test(c)) { return 'error'; }
    if (/shiny-notification-warning/.test(c)) { return 'warning'; }
    if (/shiny-notification-message/.test(c)) { return 'message'; }
    return 'default';
  }
  function scan() {
    var nodes = document.querySelectorAll('.shiny-notification');
    for (var i = 0; i < nodes.length; i++) {
      var n = nodes[i];
      var t = txt(n.querySelector('.shiny-notification-content-text') || n);
      if (!t || n.__k233 === t) { continue; }
      n.__k233 = t;
      K.notes.push({ t: Date.now(), text: t.slice(0, 4000), type: kind(n.className || ''), id: n.id || '' });
    }
    var m = document.getElementById('shiny-modal');
    var shown = !!(m && window.getComputedStyle(m).display !== 'none');
    var cur = K.modals.length ? K.modals[K.modals.length - 1] : null;
    if (shown) {
      var title = txt(m.querySelector('.modal-title'));
      var progress = !!m.querySelector('.progress, .fa-spin, .spinner, .spinner-border, .spinner-grow');
      if (!cur || cur.end !== null || cur.title !== title) {
        if (cur && cur.end === null) { cur.end = Date.now(); }
        K.modals.push({ start: Date.now(), end: null, title: title, progress: progress, text: txt(m).slice(0, 2000) });
      } else if (progress) { cur.progress = true; }
    } else if (cur && cur.end === null) { cur.end = Date.now(); }
  }
  K.timer = setInterval(function () { try { scan(); } catch (e) {} }, 200);
  scan();
  return true;
})())---"

K233_STATE_JS <- r"---((function () {
  var K = window.__k233 || { notes: [], modals: [] };
  var ids = %s, sels = %s, texts = {}, vis = [];
  ids.forEach(function (i) { var e = document.getElementById(i); texts[i] = e ? (e.innerText || e.textContent || '').trim() : ''; });
  sels.forEach(function (s) {
    var e = document.querySelector(s);
    vis.push(!!e && (typeof e.checkVisibility === 'function' ? e.checkVisibility({ checkOpacity: true, checkVisibilityCSS: true }) : e.offsetParent !== null));
  });
  var dg = document.getElementById('test_state_digest');
  var steps = Array.prototype.map.call(document.querySelectorAll('#proteomics_workflow-workflow_progress .stepper-step'), function (s) { return s.className; });
  var conn = !!(window.Shiny && Shiny.shinyapp && Shiny.shinyapp.isConnected && Shiny.shinyapp.isConnected());
  return { now: Date.now(), conn: conn, overlay: !!document.getElementById('shiny-disconnected-overlay'),
    notes: K.notes.filter(function (n) { return n.t >= %.0f; }), modals: K.modals,
    digest: dg ? dg.textContent : '', steps: steps, texts: texts, vis: vis };
})())---"

K233_TAB_JS <- r"---((function (tabset, value, click) {
  var root = tabset.charAt(0) === '[' ? document.querySelector(tabset) : document.getElementById(tabset);
  if (!root) { return { found: false, reason: 'no tabset ' + tabset }; }
  var norm = function (s) { return (s || '').replace(/\s+/g, ' ').trim(); };
  var anchors = Array.prototype.slice.call(root.querySelectorAll('a[data-toggle="tab"], a[data-bs-toggle="tab"], a[data-value]'));
  var a = anchors.filter(function (x) { return x.getAttribute('data-value') === value; })[0] ||
    anchors.filter(function (x) { return norm(x.textContent) === value; })[0];
  if (!a) { return { found: false, reason: 'no tab ' + value, titles: anchors.map(function (x) { return norm(x.textContent); }) }; }
  var on = function (x) { return x.classList.contains('active') || (!!x.parentElement && x.parentElement.classList.contains('active')); };
  var clicked = false;
  if (click && !on(a)) { a.click(); clicked = true; }
  var act = anchors.filter(on)[0];
  return { found: true, title: norm(a.textContent), active: act ? norm(act.textContent) : '', ok: !!act && act === a, clicked: clicked };
})(%s, %s, %s))---"

K233_TITLES_JS <- r"---((function (tabset) {
  var root = document.getElementById(tabset);
  if (!root) { return []; }
  return Array.prototype.map.call(root.querySelectorAll('a[data-toggle="tab"], a[data-bs-toggle="tab"]'), function (x) {
    return (x.textContent || '').replace(/\s+/g, ' ').trim(); });
})(%s))---"

K233_PRESS_JS <- r"---((function (sel, click) {
  var el = document.querySelector(sel);
  if (!el) { return { found: false }; }
  var test = 'checkVisibility', vis;
  if (typeof el.checkVisibility === 'function') { vis = el.checkVisibility({ checkOpacity: true, checkVisibilityCSS: true }); }
  else { test = 'offsetParent'; vis = el.offsetParent !== null; }
  var dis = !!el.disabled || el.classList.contains('disabled') || el.getAttribute('aria-disabled') === 'true';
  if (!vis || dis || !click) { return { found: true, visible: vis, disabled: dis, test: test }; }
  el.scrollIntoView({ block: 'center' });
  var t = Date.now();
  el.click();
  return { found: true, visible: vis, disabled: dis, test: test, clicked: true, t: t };
})(%s, %s))---"

K233_MARK_PANE_JS <- r"---((function () {
  document.querySelectorAll('[data-k233-shot]').forEach(function (e) { e.removeAttribute('data-k233-shot'); });
  var ul = document.getElementById('proteomics_workflow-workflow_tabs');
  var a = ul ? ul.querySelector('li.active > a, a.active') : null;
  var href = a ? (a.getAttribute('href') || a.getAttribute('data-bs-target')) : null;
  var p = href ? document.querySelector(href) : null;
  if (!p) { return null; }
  p.setAttribute('data-k233-shot', '1');
  var r = p.getBoundingClientRect();
  return { width: r.width, height: r.height, scroll_height: p.scrollHeight, dpr: window.devicePixelRatio };
})())---"

K233_OUTPUTS_JS <- r"---((function (ids) {
  var pane = document.querySelector('[data-k233-shot="1"]');
  var pr = pane ? pane.getBoundingClientRect() : { top: 0, left: 0 };
  return ids.map(function (i) {
    var e = document.getElementById(i);
    if (!e) { return { id: i, exists: false, kind: 'missing', nonempty: false, text: '' }; }
    var k = 'text', ne = false, text = (e.innerText || e.textContent || '').trim();
    if (e.classList.contains('shiny-plot-output') || e.classList.contains('shiny-image-output')) {
      k = 'plot'; var im = e.querySelector('img'); ne = !!(im && im.naturalWidth > 0);
    } else if (e.classList.contains('datatables') || e.querySelector('table')) {
      k = 'table'; var rows = e.querySelectorAll('tbody tr'); ne = rows.length > 0 && !e.querySelector('td.dataTables_empty');
      text = Array.prototype.slice.call(e.querySelectorAll('thead tr, tbody tr'), 0, 7).map(function (r) { return (r.innerText || '').trim(); }).join('\n');
    } else { ne = text.length > 0; }
    var r = e.getBoundingClientRect();
    return { id: i, exists: true, kind: k, nonempty: ne, recalculating: e.classList.contains('recalculating'),
      text: text.slice(0, 5000), top: r.top - pr.top, bottom: r.bottom - pr.top, left: r.left - pr.left, right: r.right - pr.left };
  });
})(%s))---"

k233Js <- function(S, script) {
    if (is.null(S$driver)) return(NULL)
    S$js_error <- NULL
    tryCatch(S$driver$get_js(script, timeout = 20000), error = function(e) {
        S$js_error <- conditionMessage(e)
        NULL
    })
}

k233Json <- function(x) jsonlite::toJSON(x, auto_unbox = TRUE)

k233WaitFor <- function(S, script, seconds) {
    started <- Sys.time()
    repeat {
        if (isTRUE(k233Js(S, script))) return(TRUE)
        if (k233Seconds(started) > seconds || k233Remaining(S$G) <= 0) return(FALSE)
        Sys.sleep(0.5)
    }
}

k233Idle <- function(S, seconds) {
    if (!is.null(S$driver)) tryCatch(S$driver$wait_for_idle(duration = 200, timeout = seconds * 1000), error = \(e) NULL)
    invisible(NULL)
}

k233State <- function(S, since, ids = character(), sels = character()) {
    script <- sprintf(K233_STATE_JS, jsonlite::toJSON(as.character(ids)), jsonlite::toJSON(as.character(sels)),
        as.numeric(since))
    st <- k233Js(S, script)
    if (is.null(st)) return(NULL)
    st$dg <- if (is.character(st$digest) && nzchar(st$digest)) {
        tryCatch(jsonlite::fromJSON(st$digest, simplifyVector = FALSE), error = \(e) NULL)
    }
    st
}

k233DigestStatus <- function(st, key) {
    value <- st$dg$step_status_per_omic$proteomics[[key]]
    if (is.null(value)) NA_character_ else as.character(value)
}

k233DigestState <- function(st) {
    value <- st$dg$r6_current_state_per_omic$proteomics
    if (is.null(value)) NA_character_ else as.character(value)
}

k233DigestHas <- function(st, state) {
    isTRUE(k233DigestState(st) == state) || state %in% as.character(unlist(st$dg$r6_state_history_per_omic$proteomics))
}

k233Stepper <- function(st, i) length(st$steps) >= i && grepl("step-completed", st$steps[[i]], fixed = TRUE)

k233Noted <- function(st, pattern, type = NULL) {
    any(vapply(st$notes, \(n) grepl(pattern, n$text, ignore.case = TRUE) && (is.null(type) || identical(n$type, type)),
        logical(1)))
}

k233Text <- function(st, id) {
    value <- st$texts[[id]]
    if (is.null(value)) "" else as.character(value)
}

k233InputValue <- function(S, id) {
    k233Js(S, sprintf(paste0("(function(id){var e=document.getElementById(id);if(!e){return null;}",
        "if(['SELECT','INPUT','TEXTAREA'].indexOf(e.tagName)>=0){return e.type==='checkbox'?e.checked:e.value;}",
        "var c=e.querySelector('input:checked');return c?c.value:null;})(%s)"), k233Json(id)))
}

k233LogAction <- function(S, step, input, value) {
    G <- S$G
    G$actions[[length(G$actions) + 1L]] <- list(session = S$label, step = step, input = input,
        value = k233Clean(value), utc = k233Utc())
}

k233Set <- function(S, step, id, value) {
    k233LogAction(S, step, id, value)
    tryCatch({
        S$G$helpers$e2e_set_input_and_idle(S$driver, id, value, timeout = 15000)
        "set"
    }, error = \(e) paste("set with an idle wait that did not settle:", conditionMessage(e)))
}

k233Selectize <- function(S, step, id, values, multiple = TRUE) {
    k233LogAction(S, step, id, paste(values, collapse = ","))
    tryCatch({
        S$G$helpers$e2e_set_selectize_values(S$driver, id, values, timeout = 30000, multiple = multiple)
        TRUE
    }, error = \(e) FALSE)
}

k233Options <- function(S, id) {
    tryCatch(S$G$helpers$e2e_selectize_option_values(S$driver, id, timeout = 30000), error = \(e) character())
}

#' Make a tab active the way a user does (a click on its link, only if it is not active) and report whether the link
#' that is active is the requested one, by the title read from the rendered tabset (5.3 (0)).
k233ActivateTab <- function(S, tabset, value) {
    k233Js(S, sprintf(K233_TAB_JS, k233Json(tabset), k233Json(value), "true"))
    started <- Sys.time()
    repeat {
        r <- k233Js(S, sprintf(K233_TAB_JS, k233Json(tabset), k233Json(value), "false"))
        if (isTRUE(r$ok)) return(sprintf("%s = %s: active", tabset, r$title))
        if (k233Seconds(started) > 15) break
        Sys.sleep(0.5)
    }
    if (!isTRUE(r$found)) return(sprintf("%s = %s: NOT FOUND (%s; titles: %s)", tabset, value, r$reason %||% "",
        paste(unlist(r$titles), collapse = ", ")))
    sprintf("%s = %s: NOT active (active: %s)", tabset, r$title %||% value, r$active %||% "")
}

k233Selector <- function(S, target, by) {
    switch(by, testid = S$G$helpers$e2e_selector(target), css = target, S$G$helpers$e2e_input_selector(target))
}

#' One press: the target must be visible and enabled (checkVisibility with opacity and visibility CSS, offsetParent
#' only where checkVisibility is missing), waited for at most 30 s while the page renders; then one click.
k233Press <- function(S, step, target, by = "id") {
    selector <- k233Selector(S, target, by)
    started <- Sys.time()
    answered <- 0
    last <- NULL
    silent <- 0L
    repeat {
        polled <- Sys.time()
        r <- k233Js(S, sprintf(K233_PRESS_JS, k233Json(selector), "false"))
        if (is.null(r)) {
            silent <- silent + 1L
        } else {
            last <- r
            answered <- answered + k233Seconds(polled) + 0.5
            if (isTRUE(r$visible) && !isTRUE(r$disabled)) break
        }
        if (answered > 30 || k233Seconds(started) > S$hard_s || k233Remaining(S$G) <= 0) {
            reason <- if (is.null(last)) "not answered by the page" else if (!isTRUE(last$found)) "not in the page" else
                if (!isTRUE(last$visible)) "not visible" else "disabled"
            k233KeepPane(S, step)
            return(list(ok = FALSE, test = last$test, reason = sprintf("%s is %s (%s; %d unanswered polls; last script error: %s)",
                target, reason, last$test %||% "", silent, S$js_error %||% "none")))
        }
        Sys.sleep(0.5)
    }
    r <- k233Js(S, sprintf(K233_PRESS_JS, k233Json(selector), "true"))
    k233LogAction(S, step, target, "click")
    if (!isTRUE(r$clicked)) return(list(ok = FALSE, reason = paste(target, "could not be clicked"), test = r$test))
    list(ok = TRUE, t = r$t, test = r$test)
}

k233KeepPane <- function(S, step) {
    html <- k233Js(S, paste0("(function(){var u=document.getElementById('proteomics_workflow-workflow_tabs');",
        "var a=u?u.querySelector('li.active > a, a.active'):null;var h=a?(a.getAttribute('href')||''):'';",
        "var p=h?document.querySelector(h):document.body;return p?p.outerHTML:'';})()"))
    writeLines(as.character(html %||% paste("no answer:", S$js_error %||% "")),
        file.path(S$G$out, sprintf("%s-%s-unreachable-dom.txt", S$label, step)))
}

k233Upload <- function(S, step, id, path, timeout_s = S$hard_s) {
    k233LogAction(S, step, id, paste("upload", basename(path)))
    error <- tryCatch({
        S$G$helpers$e2e_upload_file(S$driver, id, path, wait = FALSE, timeout = 60000)
        ""
    }, error = \(e) conditionMessage(e))
    script <- sprintf(paste0("(function(){var p=document.querySelector(%s);",
        "return p?(p.innerText||p.textContent||'').trim():'';})()"), k233Json(sprintf("#%s_progress .progress-bar", id)))
    started <- Sys.time()
    repeat {
        text <- as.character(k233Js(S, script) %||% "")
        if (grepl("complete", text, ignore.case = TRUE)) {
            return(list(ok = TRUE, text = text, seconds = round(k233Seconds(started), 1), error = error))
        }
        if (grepl("exceed|error|fail", text, ignore.case = TRUE)) {
            return(list(ok = FALSE, text = text, seconds = round(k233Seconds(started), 1), error = error))
        }
        if (k233Seconds(started) > timeout_s || k233Remaining(S$G) <= 0) {
            return(list(ok = FALSE, text = paste("upload not complete:", text, error), error = error,
                seconds = round(k233Seconds(started), 1)))
        }
        Sys.sleep(1)
    }
}

# --- one step: press, wait, classify, capture, record (5.3) ------------------------------------------------------

#' The classifier (5.3 (iii)): the first of the completion signal; a terminal signal (an error notification, a refusal
#' text in the status output, or a modal that is neither expected nor a progress modal) followed by no completion
#' within 60 s; the session's disconnection or the application's exit; the hard wait. Warnings never end a step;
#' progress modals are recorded and never end a step.
k233Wait <- function(S, spec, since) {
    started <- Sys.time()
    hard <- spec$hard_s %||% S$hard_s
    err <- NULL
    ref <- NULL
    mod <- NULL
    handled <- numeric()
    fails <- 0L
    dom_t <- NA_real_
    dig_t <- NA_real_
    last <- NULL
    result <- function(outcome, message = "") {
        list(outcome = outcome, message = message, dom_t = dom_t, dig_t = dig_t, st = last, since = since)
    }
    repeat {
        el <- k233Seconds(started)
        if (k233Remaining(S$G) <= 0) {
            return(result("deadline", "the driver's internal deadline was reached while this step was running"))
        }
        if (is.function(spec$hook)) spec$hook(S, el)
        if (!S$proc$is_alive()) {
            return(result("crashed", sprintf("the application process exited (status %s)",
                k233Clean(tryCatch(S$proc$get_exit_status(), error = \(e) "?")))))
        }
        st <- k233State(S, since, unique(c(spec$texts, spec$status_out)), spec$sels %||% character())
        if (is.null(st)) {
            fails <- fails + 1L
            if (fails >= 10L) return(result("crashed", paste("the browser page stopped answering:", S$js_error %||% "")))
            Sys.sleep(S$poll)
            next
        }
        fails <- 0L
        last <- st
        if (!isTRUE(st$conn) || isTRUE(st$overlay)) {
            return(result("crashed", "the browser session was disconnected from the application"))
        }
        open <- Filter(\(m) is.null(m$end), st$modals)
        open <- if (length(open)) open[[length(open)]] else NULL
        expected <- !is.null(open) && open$title %in% spec$modals
        if (expected && !any(handled == open$start) && is.function(spec$on_modal)) {
            handled <- c(handled, open$start)
            spec$on_modal(S, open)
        }
        d <- isTRUE(spec$dom(st))
        g <- if (S$test_mode && is.function(spec$digest)) isTRUE(spec$digest(st)) else NA
        if (d && is.na(dom_t)) dom_t <- el
        if (isTRUE(g) && is.na(dig_t)) dig_t <- el
        # A user answers an expected dialog before going on: an open expected dialog holds the step's completion.
        if (!expected && (if (is.na(g)) d else g)) return(k233Settle(S, spec, result("complete"), started))
        errors <- Filter(\(n) identical(n$type, "error"), st$notes)
        if (length(errors) && is.null(err)) err <- list(t = el, msg = errors[[1L]]$text)
        if (!is.null(spec$status_out) && !is.null(spec$refusal) && is.null(ref)) {
            text <- k233Text(st, spec$status_out)
            if (grepl(spec$refusal, text, ignore.case = TRUE)) ref <- list(t = el, msg = text)
        }
        if (is.null(open)) {
            mod <- NULL
        } else if (expected) {
            mod <- NULL
        } else if (!isTRUE(open$progress) && (is.null(mod) || !identical(mod$start, open$start))) {
            mod <- list(t = el, start = open$start, msg = sprintf("modal '%s': %s", open$title, open$text))
        }
        for (signal in list(err, ref, mod)) {
            if (!is.null(signal) && el - signal$t >= 60) return(result("refused", signal$msg))
        }
        if (el >= hard) return(result("timeout", sprintf("no completion signal within %d s", as.integer(hard))))
        Sys.sleep(S$poll)
    }
}

#' After a test-mode completion read from the digest, the DOM signal is given up to 30 s (L-5 compares the two), then
#' the page is let settle, as a user waits for the busy indicator to clear.
k233Settle <- function(S, spec, res, started) {
    if (S$test_mode && is.function(spec$digest) && is.na(res$dom_t)) {
        limit <- k233Seconds(started) + 30
        repeat {
            el <- k233Seconds(started)
            if (el > limit || k233Remaining(S$G) <= 0) break
            st <- k233State(S, res$since, unique(c(spec$texts, spec$status_out)), spec$sels %||% character())
            if (!is.null(st)) {
                res$st <- st
                if (isTRUE(spec$dom(st))) {
                    res$dom_t <- el
                    break
                }
            }
            Sys.sleep(S$poll)
        }
    }
    res$agreement <- if (!S$test_mode || !is.function(spec$digest)) {
        "dom_only"
    } else if (!is.na(res$dom_t) && !is.na(res$dig_t)) {
        "agree"
    } else {
        "disagree"
    }
    k233Idle(S, 15)
    res
}

#' One action of a step (5.3): tabs, settings, the reachability check, one click, the wait, the classification, the
#' captures and one row of steps.tsv; findings follow from the row. `spec$custom` replaces the press and the wait
#' (the report download).
k233Action <- function(S, spec) {
    G <- S$G
    if (k233Remaining(G) <= 0) return(k233NotReachedRow(S, spec$step, spec$name, "not_reached (deadline)"))
    S$nn <- S$nn + 1L
    spec$hook <- spec$hook %||% S$hooks[[spec$step]]
    started <- Sys.time()
    rlog0 <- length(k233ReadLog(S$log))
    since <- as.numeric(k233Js(S, "Date.now()") %||% 0)
    checks <- vapply(spec$tabs %||% list(), \(tb) k233ActivateTab(S, tb[[1L]], tb[[2L]]), character(1))
    prepared <- if (is.function(spec$prepare)) {
        tryCatch(spec$prepare(S), error = \(e) list(refuse = paste("driver error while setting inputs:",
            conditionMessage(e))))
    }
    press <- list(test = "")
    if (is.list(prepared) && !is.null(prepared$refuse)) {
        res <- list(outcome = "refused", message = prepared$refuse, unreachable = TRUE, since = since)
    } else {
        if (length(spec$tabs)) checks <- vapply(spec$tabs, \(tb) k233ActivateTab(S, tb[[1L]], tb[[2L]]), character(1))
        if (is.function(spec$custom)) {
            res <- spec$custom(S)
            res$since <- since
        } else {
            press <- k233Press(S, spec$step, spec$target, spec$by %||% "id")
            res <- if (isTRUE(press$ok)) {
                k233Wait(S, spec, press$t)
            } else {
                list(outcome = "refused", message = paste("not reachable as a user:", press$reason), unreachable = TRUE,
                    since = since)
            }
        }
    }
    seconds <- round(k233Seconds(started), 1)
    cap <- k233Capture(S, spec, res, since)
    modals <- Filter(\(m) as.numeric(m$start) >= since, res$st$modals %||% cap$modals)
    shown <- function(m) sprintf("%s (%.0f s)", m$title,
        ((if (is.null(m$end)) as.numeric(res$st$now %||% m$start) else as.numeric(m$end)) - as.numeric(m$start)) / 1000)
    warnings <- vapply(Filter(\(n) identical(n$type, "warning"), cap$notes), \(n) n$text, character(1))
    row <- list(
        session = S$label, step = spec$step, name = spec$name, outcome = res$outcome, seconds = seconds,
        started_utc = k233Utc(started), message = res$message %||% "",
        head_state = k233DigestState(res$st %||% list()),
        tab_checks = if (length(checks)) paste(checks, collapse = "; ") else spec$tab_note %||% "none (no tab)",
        visibility_test = press$test %||% "", digest_signal_s = res$dig_t %||% NA, dom_signal_s = res$dom_t %||% NA,
        agreement = if (identical(res$outcome, "complete")) res$agreement %||% "dom_only" else "n/a",
        progress_modals = paste(vapply(Filter(\(m) isTRUE(m$progress), modals), shown, character(1)), collapse = "; "),
        other_modals = paste(vapply(Filter(\(m) !isTRUE(m$progress), modals), shown, character(1)), collapse = "; "),
        warnings = paste(warnings, collapse = " || "), notifications = length(cap$notes),
        console_errors = length(cap$console_errors), named_outputs = cap$outputs_summary,
        files = paste(sprintf("%s %s", names(cap$files), cap$files), collapse = " | "),
        rlog_lines = sprintf("%s:%d-%d", basename(S$log), rlog0 + 1L, length(k233ReadLog(S$log))),
        note = paste(c(res$note, cap$note), collapse = "; ")
    )
    G$steps[[length(G$steps) + 1L]] <- row
    k233StepFindings(S, spec, res, cap, row)
    if (identical(res$outcome, "complete") && is.function(spec$after)) {
        tryCatch(spec$after(S, res, cap), error = \(e) k233Observation(G, S$label, paste(spec$step, "after"),
            paste("driver error:", conditionMessage(e))))
    }
    k233WriteAll(G)
    res$row <- row
    res$cap <- cap
    res
}

k233NotReachedRow <- function(S, step, name, outcome = "not_reached") {
    S$G$steps[[length(S$G$steps) + 1L]] <- list(session = S$label, step = step, name = name, outcome = outcome,
        seconds = 0, started_utc = k233Utc(), message = "", tab_checks = "", agreement = "n/a")
    list(outcome = outcome)
}

k233CheckOutputs <- function(S, ids) {
    if (!length(ids)) return(list())
    k233Js(S, sprintf(K233_OUTPUTS_JS, jsonlite::toJSON(as.character(ids)))) %||% list()
}

#' Captures after each action (5.3 (v)): the active tab pane's element screenshot (and per-output element screenshots
#' where the pane image does not cover an output), a viewport screenshot, the notifications, the console, the named
#' outputs' text and state, and the digest's head and history in test mode.
k233Capture <- function(S, spec, res, since) {
    G <- S$G
    stem <- sprintf("%s-%02d-%s", S$label, S$nn, spec$step)
    files <- character()
    keep <- function(path) {
        if (file.exists(path) && file.size(path) > 0) files[[basename(path)]] <<- k233Sha(path)
    }
    note <- character()
    if (is.null(S$driver)) return(list(files = files, notes = list(), console_errors = character(), modals = list(),
        outputs_summary = "", outputs = list(), note = "no browser session"))
    box <- k233Js(S, K233_MARK_PANE_JS)
    pane <- file.path(G$out, paste0(stem, ".png"))
    shot <- tryCatch({
        S$driver$get_screenshot(pane, selector = if (is.list(box)) "[data-k233-shot='1']" else "body")
        "ok"
    }, error = \(e) conditionMessage(e))
    if (!identical(shot, "ok")) note <- c(note, paste("pane screenshot:", shot))
    keep(pane)
    view <- file.path(G$out, paste0(stem, "-viewport.png"))
    tryCatch(S$driver$get_screenshot(view, selector = "viewport"), error = \(e) NULL)
    keep(view)
    outputs <- k233CheckOutputs(S, spec$outputs)
    size <- k233PngSize(pane)
    ratio <- if (is.list(box) && is.finite(size[[1L]]) && box$width > 0) size[[1L]] / box$width else NA_real_
    for (o in outputs) {
        covered <- isTRUE(o$exists) && is.finite(ratio) && o$bottom * ratio <= size[[2L]] + 2 &&
            o$right * ratio <= size[[1L]] + 2 && o$top >= -1
        zero <- isTRUE(o$exists) && isTRUE(o$bottom - o$top <= 0)
        how <- if (covered) "pane image" else if (zero) "zero size (nothing drawn)" else "not covered"
        if (!covered && !zero && isTRUE(o$exists)) {
            element <- file.path(G$out, sprintf("%s-out-%s.png", stem, sub(".*-", "", o$id)))
            k233Js(S, sprintf("(function(){var e=document.getElementById(%s);if(e){e.scrollIntoView();}return !!e;})()",
                k233Json(o$id)))
            tryCatch(S$driver$get_screenshot(element, selector = paste0("#", o$id)), error = \(e) NULL)
            keep(element)
            if (file.exists(element) && file.size(element) > 0) how <- "own element image"
        }
        G$coverage[[length(G$coverage) + 1L]] <- list(session = S$label, step = spec$step, output = o$id,
            pane_height = if (is.list(box)) round(box$height) else NA, image = paste(size, collapse = "x"),
            covered = how)
    }
    notes <- k233Js(S, sprintf("(window.__k233 ? window.__k233.notes : []).filter(function(n){return n.t >= %.0f;})",
        since)) %||% list()
    modals <- k233Js(S, "(window.__k233 ? window.__k233.modals : [])") %||% list()
    visible <- tryCatch(G$helpers$e2e_collect_notifications(S$driver), error = \(e) list())
    notes_path <- file.path(G$out, paste0(stem, "-notifications.json"))
    jsonlite::write_json(list(during_step = notes, visible_at_end = visible, modals = modals), notes_path,
        auto_unbox = TRUE, pretty = TRUE, null = "null")
    keep(notes_path)
    logs <- tryCatch(as.data.frame(S$driver$get_logs()), error = \(e) data.frame())
    fresh <- if (nrow(logs) > S$console_seen) logs[seq.int(S$console_seen + 1L, nrow(logs)), , drop = FALSE] else logs[0, ]
    S$console_seen <- nrow(logs)
    console_path <- file.path(G$out, paste0(stem, "-console.json"))
    jsonlite::write_json(fresh, console_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
    keep(console_path)
    level <- tolower(as.character(fresh$level %||% rep("", nrow(fresh))))
    message <- as.character(fresh$message %||% rep("", nrow(fresh)))
    location <- as.character(fresh$location %||% rep("chromote", nrow(fresh)))
    bad <- (level %in% c("error", "severe") |
        grepl("\\b(uncaught|typeerror|referenceerror|syntaxerror|rangeerror|securityerror)\\b", message,
            ignore.case = TRUE)) & !grepl("shiny", location, ignore.case = TRUE)
    st <- res$st
    outputs_path <- file.path(G$out, paste0(stem, "-outputs.json"))
    jsonlite::write_json(list(step = spec$step, outcome = res$outcome, named_outputs = outputs,
        status_texts = st$texts, stepper = st$steps,
        digest_head = if (S$test_mode) k233DigestState(st %||% list()) else NULL,
        digest_history = if (S$test_mode) st$dg$r6_state_history_per_omic$proteomics else NULL,
        digest_tab_status = if (S$test_mode) st$dg$step_status_per_omic$proteomics else NULL),
        outputs_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
    keep(outputs_path)
    summary <- paste(vapply(outputs, \(o) sprintf("%s=%s", sub(K233_NS, "", o$id),
        if (isTRUE(o$nonempty)) "ok" else if (isTRUE(o$exists)) "EMPTY" else "MISSING"), character(1)), collapse = " ")
    list(files = files, notes = notes, modals = modals, console_errors = message[bad], outputs = outputs,
        outputs_summary = summary, note = paste(note, collapse = "; "))
}

k233Finding <- function(S, step, kind, message, row = NULL, planted = FALSE, extra = list()) {
    G <- S$G
    acts <- Filter(\(a) identical(a$session, S$label) && identical(a$step, step), G$actions)
    shots <- if (is.null(row)) "" else grep("\\.png ", strsplit(k233Clean(row$files), " \\| ")[[1L]], value = TRUE)
    rlog <- k233ReadLog(S$log)
    range <- if (is.null(row)) integer() else {
        bounds <- as.integer(strsplit(sub("^.*:", "", row$rlog_lines), "-")[[1L]])
        if (length(bounds) == 2L && bounds[[2L]] >= bounds[[1L]]) seq.int(bounds[[1L]], bounds[[2L]]) else integer()
    }
    lines <- rlog[range[range <= length(rlog)]]
    excerpt <- grep("ERROR|Error|WARN|error", lines, value = TRUE)
    excerpt <- utils::head(if (length(excerpt)) excerpt else utils::tail(lines, 3L), 3L)
    G$findings[[length(G$findings) + 1L]] <- c(list(
        id = sprintf("K233-F%02d", length(G$findings) + 1L), session = S$label,
        run_id = Sys.getenv("GITHUB_RUN_ID"), commit = G$phase$commit, input = G$phase$inputs$id %||% "",
        input_sha256 = G$phase$inputs$sha %||% "", step = step, kind = kind, planted = planted,
        actions = paste(vapply(acts, \(a) sprintf("%s=%s@%s", a$input, a$value, a$utc), character(1)), collapse = "; "),
        screenshot = paste(shots, collapse = "; "), message = message, rlog = row$rlog_lines %||% "",
        rlog_excerpt = paste(substr(excerpt, 1L, 300L), collapse = " || "),
        head_state = row$head_state %||% "", known = "classified in outcome.md"
    ), extra)
    invisible(NULL)
}

k233StepFindings <- function(S, spec, res, cap, row) {
    planted <- isTRUE(spec$planted)
    if (res$outcome %in% c("refused", "timeout", "crashed")) {
        k233Finding(S, spec$step, 1L, sprintf("%s: %s", res$outcome, res$message %||% ""), row, planted)
    }
    if (identical(res$outcome, "deadline")) k233Finding(S, spec$step, 5L, res$message, row, planted)
    if (row$seconds > S$G$bounds$finding_s) {
        k233Finding(S, spec$step, 5L, sprintf("step took %.0f s (bound %d s)", row$seconds, S$G$bounds$finding_s), row)
    }
    if (!identical(res$outcome, "complete")) return(invisible(NULL))
    for (o in cap$outputs) {
        if (!isTRUE(o$nonempty)) {
            k233Finding(S, spec$step, 3L, sprintf("named output %s: %s", o$id, if (isTRUE(o$exists)) "empty" else
                "missing"), row, planted)
        }
    }
    errors <- Filter(\(n) identical(n$type, "error"), cap$notes)
    for (n in errors) k233Finding(S, spec$step, 6L, paste("error notification on a completed step:", n$text), row)
    for (m in cap$console_errors) k233Finding(S, spec$step, 6L, paste("browser console error:", substr(m, 1L, 500L)), row)
    invisible(NULL)
}

# ---------------------------------------------------------------------------------------------------------------
# The route runner and the steps (5.2).
# ---------------------------------------------------------------------------------------------------------------
k233RunRoute <- function(S, route, R, upto = NULL, hooks = list()) {
    S$hooks <- hooks
    res <- list(outcome = "not_reached")
    for (i in seq_along(route$steps)) {
        step <- route$steps[[i]]
        res <- if (k233Remaining(S$G) <= 0) {
            k233NotReachedRow(S, step$id, step$name, "not_reached (deadline)")
        } else {
            match.fun(step$fn)(S, R, step, route)
        }
        if (identical(res$outcome, "complete") && is.function(R$after[[step$id]])) R$after[[step$id]](S, R)
        stop_here <- !identical(res$outcome, "complete") &&
            !(isTRUE(step$continue_on_refusal) && identical(res$outcome, "refused"))
        if (stop_here) {
            rest <- route$steps[-seq_len(i)]
            reason <- if (grepl("deadline", res$outcome)) "not_reached (deadline)" else "not_reached"
            for (r in rest) k233NotReachedRow(S, r$id, r$name, reason)
            return(res)
        }
        if (identical(step$id, upto)) return(res)
    }
    res
}

k233Workflow <- function(value) list(paste0(K233_NS, "workflow_tabs"), value)

k233ProjectPaths <- function(R) {
    base <- file.path(R$project_dir, R$label)
    list(base = base, source = file.path(base, "scripts", "proteomics"),
        da = file.path(base, "results", "proteomics", "da_proteins"))
}

k233StepHome <- function(S, R, step, route) {
    prepare <- function(S) {
        for (t in c("tile-proteomics", "btn-start-analysis")) {
            p <- k233Press(S, "S1", t, "testid")
            if (!isTRUE(p$ok)) return(list(refuse = paste("not reachable as a user:", p$reason)))
            k233Idle(S, 10)
        }
        if (!k233WaitFor(S, "document.getElementById('experiment_label') !== null", 30)) {
            return(list(refuse = "the Experiment Setup dialog did not open within 30 s"))
        }
        k233Set(S, "S1", "experiment_label", R$label)
        k233Set(S, "S1", "project_dir", R$project_dir)
        NULL
    }
    k233Action(S, list(step = "S1", name = step$name, target = "btn-confirm-setup", by = "testid", prepare = prepare,
        tab_note = "none (home page and the Experiment Setup dialog)", sels = "[data-testid='prot-tab-setup']",
        dom = \(st) isTRUE(st$vis[[1L]]), digest = \(st) "proteomics" %in% unlist(st$dg$project_dir_keys),
        modals = "Experiment Setup", outputs = paste0(K233_NS, "workflow_progress")))
}

k233ChooseFile <- function(S, step, button, path) {
    p <- k233Press(S, step, button)
    if (!isTRUE(p$ok)) return(list(ok = FALSE, how = paste("chooser button:", p$reason)))
    modal <- sprintf("#%s-modal", button)
    if (!k233WaitFor(S, sprintf("document.querySelectorAll(%s).length > 0", k233Json(paste(modal, ".sF-fileList > *"))),
        30)) {
        return(list(ok = FALSE, how = "the file chooser listed nothing within 30 s"))
    }
    shot <- file.path(S$G$out, sprintf("%s-%s-chooser-%s.png", S$label, step, sub(".*-", "", button)))
    tryCatch(S$driver$get_screenshot(shot, selector = paste(modal, ".modal-content")), error = \(e) NULL)
    parts <- c(strsplit(dirname(path), "/", fixed = TRUE)[[1L]][-1L], basename(path))
    find <- function(name) sprintf(paste0("(function(m,n){var f=Array.prototype.filter.call(document.querySelectorAll(m+",
        "' .sF-fileList .sF-file'),function(x){var d=x.querySelector('.sF-file-name>div');return d&&d.textContent===n;})",
        "[0];return f||null;})(%s,%s)"), k233Json(modal), k233Json(name))
    for (i in seq_along(parts)) {
        exists <- sprintf("(%s) !== null", find(parts[[i]]))
        if (!k233WaitFor(S, exists, 20)) return(list(ok = FALSE, how = paste("the chooser does not list", parts[[i]])))
        event <- if (i < length(parts)) "dblclick" else "click"
        k233Js(S, sprintf("(function(){var f=%s;if(!f){return false;}jQuery(f).trigger('%s');return true;})()",
            find(parts[[i]]), event))
        k233LogAction(S, step, paste(button, "chooser"), paste(event, parts[[i]]))
    }
    selector <- paste(modal, "#sF-selectButton")
    p <- k233Press(S, step, selector, "css")
    if (!isTRUE(p$ok)) return(list(ok = FALSE, how = paste("Select:", p$reason)))
    ok <- k233WaitFor(S, sprintf("document.querySelector(%s) === null", k233Json(modal)), 20)
    list(ok = ok, how = "clicked through the chooser's folder list")
}

k233StepImport <- function(S, R, step, route) {
    ns <- paste0(K233_NS, "setup_import-")
    prepare <- function(S) {
        if (S$test_mode) {
            up <- k233Upload(S, "S2", paste0(ns, "search_results_standard"), R$report)
            S$import_note <- sprintf("search results upload %s in %s s", up$text, up$seconds)
            if (!up$ok) return(list(refuse = paste("search results upload:", up$text)))
        } else {
            chosen <- k233ChooseFile(S, "S2", paste0(ns, "search_results"), R$report)
            k233Observation(S$G, S$label, "S2 search results chooser", chosen$how)
            if (!chosen$ok) return(list(refuse = paste("not reachable as a user: search results chooser:", chosen$how)))
        }
        k233WaitFor(S, sprintf("(document.getElementById(%s)||{innerText:''}).innerText.trim().length > 0",
            k233Json(paste0(ns, "format_detection"))), 60)
        detected <- as.character(k233Js(S, sprintf("(document.getElementById(%s)||{innerText:''}).innerText.trim()",
            k233Json(paste0(ns, "format_detection")))) %||% "")
        k233Observation(S$G, S$label, "S2 detected format", detected)
        if (!grepl("DIA-?NN", detected, ignore.case = TRUE)) {
            k233Finding(S, "S2", 2L, paste("detected format is not DIA-NN; DIA-NN chosen by hand. Shown:", detected))
            k233Set(S, "S2", paste0(ns, "format_override"), "diann")
        }
        if (S$test_mode) {
            up <- k233Upload(S, "S2", paste0(ns, "fasta_file_standard"), R$fasta)
            S$import_note <- paste(S$import_note, sprintf("; FASTA upload %s in %s s", up$text, up$seconds))
            if (!up$ok) return(list(refuse = paste("FASTA upload:", up$text)))
        } else {
            chosen <- k233ChooseFile(S, "S2", paste0(ns, "fasta_file"), R$fasta)
            k233Observation(S$G, S$label, "S2 FASTA chooser", chosen$how)
            if (!chosen$ok) return(list(refuse = paste("not reachable as a user: FASTA chooser:", chosen$how)))
        }
        k233Set(S, "S2", paste0(ns, "mixed_species_fasta"), TRUE)
        NULL
    }
    on_modal <- function(S, modal) {
        S$organism_modal <- modal$text
        values <- list(selected_organism = k233InputValue(S, paste0(ns, "selected_organism")),
            filter_to_organism = k233InputValue(S, paste0(ns, "filter_to_organism")))
        p <- k233Press(S, "S2", paste0(ns, "confirm_organism"))
        S$organism_answered <- isTRUE(p$ok)
        k233Observation(S$G, S$label, "S2 organism modal", sprintf("%s | rendered values: %s | answered: %s",
            modal$text, k233Clean(paste(names(values), values, sep = "=")), isTRUE(p$ok)))
    }
    k233Action(S, list(step = "S2", name = step$name, target = "prot-import-process", by = "testid",
        tabs = list(k233Workflow("setup")), prepare = prepare,
        digest = \(st) identical(k233DigestStatus(st, "setup_import"), "complete"), dom = \(st) k233Stepper(st, 1L),
        modals = "Select Primary Organism", on_modal = on_modal,
        outputs = paste0(ns, c("format_detection", "import_status", "data_summary"))))
}

k233MatchContrast <- function(options, pair) {
    first <- k233Norm(pair[[1L]])
    second <- k233Norm(pair[[2L]])
    for (o in options) {
        n <- k233Norm(o)
        a <- regexpr(first, n, fixed = TRUE)
        b <- regexpr(second, n, fixed = TRUE)
        if (a > 0 && b > 0 && a < b) return(o)
    }
    NA_character_
}

k233StepDesign <- function(S, R, step, route) {
    builder <- function(x) paste0(K233_NS, "design_matrix-builder-", x)
    tabs <- builder("main_tabs")
    prepare <- function(S) {
        if (!k233WaitFor(S, sprintf("document.getElementById(%s) !== null", k233Json(builder("selected_runs"))), 60)) {
            return(list(refuse = "the design builder was not shown within 60 s"))
        }
        k233ActivateTab(S, tabs, "Factors")
        for (g in R$groups) {
            k233Set(S, "S3", builder("new_factor"), g)
            p <- k233Press(S, "S3", builder("add_factor"))
            if (!isTRUE(p$ok)) return(list(refuse = paste("not reachable as a user:", p$reason)))
            k233Idle(S, 10)
        }
        k233ActivateTab(S, tabs, "Assign Metadata")
        options <- k233Options(S, builder("selected_runs"))
        for (g in R$groups) {
            runs <- options[k233Norm(options) %in% k233Norm(R$members[[g]])]
            if (length(runs) != length(R$members[[g]])) {
                return(list(refuse = sprintf("driver: %d of %d runs of %s found among the builder's %d runs", length(runs),
                    length(R$members[[g]]), g, length(options))))
            }
            k233Selectize(S, "S3", builder("selected_runs"), runs)
            k233WaitFor(S, sprintf("document.getElementById(%s) !== null", k233Json(builder("replicate_start"))), 30)
            tryCatch(S$G$helpers$e2e_wait_for_select_options(S$driver, builder("factor1_select"), g, timeout = 30000),
                error = \(e) NULL)
            k233Selectize(S, "S3", builder("factor1_select"), g, multiple = FALSE)
            p <- k233Press(S, "S3", builder("assign_metadata"))
            if (!isTRUE(p$ok)) return(list(refuse = paste("not reachable as a user:", p$reason)))
            k233Idle(S, 10)
        }
        k233ActivateTab(S, tabs, "Contrasts")
        for (pair in R$contrasts) {
            choices <- k233Options(S, builder("contrast_group1"))
            matched <- S$G$helpers$e2e_match_available_values(pair, choices)
            k233Selectize(S, "S3", builder("contrast_group1"), matched[[1L]], multiple = FALSE)
            k233Selectize(S, "S3", builder("contrast_group2"), matched[[2L]], multiple = FALSE)
            p <- k233Press(S, "S3", builder("add_contrast"))
            if (!isTRUE(p$ok)) return(list(refuse = paste("not reachable as a user:", p$reason)))
            k233Idle(S, 10)
        }
        NULL
    }
    k233Action(S, list(step = "S3", name = step$name, target = builder("save_results"),
        tabs = list(k233Workflow("design")), prepare = prepare,
        digest = \(st) identical(k233DigestStatus(st, "design_matrix"), "complete"), dom = \(st) k233Stepper(st, 2L),
        outputs = paste0(K233_NS, "design_matrix-", c("design_matrix_preview", "contrasts_preview"))))
}

k233QcTabs <- function(route, group) {
    live <- tryCatch(MultiScholaR:::omicGuiConfig("proteomics")$qc$groups, error = \(e) NULL)
    declared <- route$qc[[group]]
    for (g in live) {
        if (!identical(g$id, group)) next
        tabs <- Filter(\(t) route$id %in% unlist(t$routes), g$tabs)
        titles <- vapply(tabs, \(t) t$title, character(1))
        ids <- vapply(tabs, \(t) t$id, character(1))
        return(list(source = "live omicGuiConfig", map = stats::setNames(ids, titles), declared = declared))
    }
    list(source = "route table (live config not read)",
        map = stats::setNames(vapply(declared, \(t) t$id, character(1)), vapply(declared, \(t) t$title, character(1))),
        declared = declared)
}

k233StepQc <- function(S, R, step, route) {
    qc <- function(x) paste0(K233_NS, "quality_control-", x)
    group_tabs <- list(k233Workflow("qc"), list(qc(route$qc_tabset), step$title))
    for (tb in group_tabs) k233ActivateTab(S, tb[[1L]], tb[[2L]])
    titles <- character()
    started <- Sys.time()
    while (!length(titles) && k233Seconds(started) < 60) {
        titles <- as.character(unlist(k233Js(S, sprintf(K233_TITLES_JS, k233Json(qc(step$tabset))))))
        if (!length(titles)) Sys.sleep(1)
    }
    tabs <- k233QcTabs(route, step$group)
    k233Observation(S$G, S$label, paste(step$id, "rendered sub-tabs"), sprintf("%s (titles from %s; expected %s)",
        paste(titles, collapse = " > "), tabs$source, paste(names(tabs$map), collapse = " > ")))
    if (!length(titles)) {
        return(k233Action(S, list(step = paste0(step$id, ".1"), name = step$name, prepare = \(S) list(
            refuse = paste("not reachable as a user: the", step$title, "sub-tabs were not rendered within 60 s")))))
    }
    res <- list(outcome = "complete")
    for (i in seq_along(titles)) {
        id <- unname(tabs$map[titles[[i]]])
        spec_tab <- Filter(\(t) identical(t$id, id), tabs$declared)
        if (is.na(id) || !length(spec_tab)) {
            return(k233Action(S, list(step = sprintf("%s.%d", step$id, i), name = titles[[i]], prepare = \(S) list(
                refuse = paste("driver: rendered sub-tab not in the route table:", titles[[i]])))))
        }
        spec_tab <- spec_tab[[1L]]
        prefix <- qc(paste0(step$group, "-", id, "-"))
        state <- spec_tab$state
        spec <- list(step = sprintf("%s.%d", step$id, i), name = titles[[i]], target = paste0(prefix, spec_tab$button),
            tabs = c(group_tabs, list(list(qc(step$tabset), titles[[i]]))),
            digest = \(st) k233DigestHas(st, state), dom = \(st) k233Noted(st, "applied successfully", "message"),
            outputs = paste0(prefix, spec_tab$outputs))
        res <- k233Action(S, spec)
        if (identical(res$outcome, "refused") && isTRUE(R$permissive) && length(spec_tab$permissive)) {
            retry <- spec
            retry$step <- paste0(spec$step, "r")
            retry$name <- paste(titles[[i]], "(the lane's permissive value, L-1's rule)")
            permissive <- spec_tab$permissive
            retry$prepare <- function(S) {
                for (nm in names(permissive)) k233Set(S, retry$step, paste0(prefix, nm), permissive[[nm]])
                NULL
            }
            res <- k233Action(S, retry)
        }
        if (!identical(res$outcome, "complete")) return(res)
    }
    res
}

k233StepNormalisation <- function(S, R, step, route) {
    norm <- function(x) paste0(K233_NS, "normalization-", x)
    tabs <- list(k233Workflow("normalization"))
    prepare <- function(S) {
        S$G$ruv_panel <- k233Js(S, sprintf("document.getElementById(%s) !== null", k233Json(norm("ruv_design_summary"))))
        k233Observation(S$G, S$label, "S6 rendered defaults", sprintf("ruv_mode=%s method=%s ruv_design_panel_present=%s",
            k233Clean(k233InputValue(S, norm("ruv_mode"))), k233Clean(k233InputValue(S, norm("norm_method"))),
            k233Clean(S$G$ruv_panel)))
        NULL
    }
    first <- k233Action(S, list(step = "S6.1", name = "normalisation and RUV at the defaults", target = "prot-norm-run",
        by = "testid", tabs = tabs, prepare = prepare,
        # No digest signal: the state ruv_corrected is committed before the step ends (run L rerun 1, 45525410012).
        dom = \(st) k233Noted(st, "Normali[sz]ation (and RUV correction )?completed", "message"),
        outputs = norm(c("norm_log", "ruv_optimization_summary"))))
    if (!identical(first$outcome, "complete")) return(first)
    summary_id <- norm("correlation_filter_summary")
    k233Action(S, list(step = "S6.2", name = "correlation filter at the rendered threshold",
        target = norm("apply_correlation_filter"), tabs = c(tabs, list(list(norm("norm_qc_tabs"), "Correlation Filtering"))),
        prepare = function(S) {
            k233Observation(S$G, S$label, "S6.2 rendered threshold",
                k233Clean(k233InputValue(S, norm("min_pearson_correlation_threshold"))))
            NULL
        },
        digest = \(st) k233DigestHas(st, "correlation_filtered") &&
            identical(k233DigestStatus(st, "normalization"), "complete"),
        dom = \(st) k233Stepper(st, 4L), texts = summary_id, outputs = norm(c("correlation_filter_summary", "final_qc_plot")),
        after = function(S, res, cap) {
            text <- k233Text(res$st, summary_id)
            S$n1_shown <- list(text = text,
                proteins = suppressWarnings(as.integer(sub("(?s).*Proteins remaining: ([0-9]+).*", "\\1", text, perl = TRUE))),
                samples = suppressWarnings(as.integer(sub("(?s).*Samples remaining: ([0-9]+).*", "\\1", text, perl = TRUE))))
        }))
}

k233StepExport <- function(S, R, step, route) {
    k233Action(S, list(step = "S7", name = step$name, target = "prot-norm-export-session", by = "testid",
        tabs = list(k233Workflow("normalization")),
        dom = \(st) k233Noted(st, "Filtered session data exported successfully", "message"),
        after = function(S, res, cap) k233NumberN1(S, R)))
}

k233StepDa <- function(S, R, step, route) {
    da <- function(x) paste0(K233_NS, "differential_expression-", x)
    tabs <- list(k233Workflow("da"))
    loaded <- k233Action(S, list(step = "S8.1", name = "load the filtered session", target = "prot-da-load-session",
        by = "testid", tabs = tabs, texts = da(c("da_status", "contrasts_display")),
        dom = \(st) grepl("Session loaded successfully", k233Text(st, da("da_status")), fixed = TRUE) ||
            k233Noted(st, "loaded successfully", "message"),
        outputs = da("contrasts_display")))
    if (!identical(loaded$outcome, "complete")) return(loaded)
    run <- k233Action(S, list(step = "S8.2", name = "run DA at the defaults", target = "prot-da-run", tabs = tabs,
        prepare = function(S) {
            k233Observation(S$G, S$label, "S8 rendered defaults", sprintf("formula=%s q=%s lfc=%s",
                k233Clean(k233InputValue(S, da("formula_string"))), k233Clean(k233InputValue(S, da("da_q_val_thresh"))),
                k233Clean(k233InputValue(S, da("treat_lfc_cutoff")))))
            NULL
        },
        texts = da("da_status"), status_out = da("da_status"), refusal = "\\[ERROR\\]|^Error",
        digest = \(st) identical(k233DigestStatus(st, "differential_expression"), "complete"),
        dom = \(st) k233Stepper(st, 5L) && grepl("Analysis Complete", k233Text(st, da("da_status")), fixed = TRUE),
        outputs = da("da_status")))
    if (!identical(run$outcome, "complete")) return(run)
    k233Action(S, list(step = "S8.3", name = "open the DA results table", by = "css",
        target = sprintf("#%s a[data-value='DA Results Table']", da("da_results_tabs")), tabs = tabs,
        dom = \(st) isTRUE(st$vis[[1L]]), sels = paste0("#", da("da_results_table"), " tbody tr"),
        outputs = da(c("da_summary_stats", "da_results_table")),
        after = function(S, res, cap) k233NumberN2(S, R)))
}

k233StepEnrichment <- function(S, R, step, route) {
    enr <- function(x) paste0(K233_NS, "enrichment_analysis-", x)
    k233Action(S, list(step = "S9", name = step$name, target = enr("run_enrichment_analysis"),
        tabs = list(k233Workflow("enrichment"), list(enr("enrichment_method_tabs"), "gprofiler2")),
        prepare = function(S) {
            tryCatch(S$G$helpers$e2e_wait_for_js(S$driver, sprintf(
                "(function(){var e=document.getElementById(%s);return !!(e && e.options && e.options.length > 0);})()",
                k233Json(enr("selected_contrast"))), timeout = 60000), error = \(e) NULL)
            options <- as.character(unlist(k233Js(S, sprintf(
                "Array.prototype.map.call((document.getElementById(%s)||{options:[]}).options,function(o){return o.value;})",
                k233Json(enr("selected_contrast"))))))
            chosen <- k233MatchContrast(options, R$contrasts[[1L]])
            if (is.na(chosen) && length(options)) chosen <- options[[1L]]
            if (!is.na(chosen)) k233Set(S, "S9", enr("selected_contrast"), chosen)
            k233Observation(S$G, S$label, "S9 rendered settings", sprintf("contrast=%s taxon=%s up=%s down=%s q=%s",
                k233Clean(chosen), k233Clean(k233InputValue(S, enr("organism_taxid"))),
                k233Clean(k233InputValue(S, enr("up_cutoff"))), k233Clean(k233InputValue(S, enr("down_cutoff"))),
                k233Clean(k233InputValue(S, enr("q_cutoff")))))
            NULL
        },
        texts = enr("enrichment_status"), status_out = enr("enrichment_status"),
        refusal = "needs protein annotation|refused|^Error",
        digest = \(st) identical(k233DigestStatus(st, "enrichment_analysis"), "complete"),
        dom = \(st) grepl("Analysis Complete", k233Text(st, enr("enrichment_status")), fixed = TRUE) &&
            grepl("Method:", k233Text(st, enr("enrichment_status")), fixed = TRUE),
        outputs = enr(c("enrichment_status", "gprofiler_results_table", "gprofiler_summary_stats"))))
}

k233StepSummary <- function(S, R, step, route) {
    sm <- function(x) paste0(K233_NS, "session_summary-", x)
    tabs <- list(k233Workflow("session_summary"))
    summary_id <- sm("session_summary")
    copy_id <- sm("copy_status")
    actions <- list(
        list(step = "S10.1", name = "Save Workflow Arguments", target = "prot-summary-save-args",
            prepare = function(S) {
                k233Set(S, "S10.1", sm("experiment_label"), R$label)
                k233Set(S, "S10.1", sm("description"), R$description)
                NULL
            },
            dom = \(st) grepl("Parameters saved [OK]", k233Text(st, summary_id), fixed = TRUE),
            outputs = summary_id, after = function(S, res, cap) k233NumberN3(S, R)),
        list(step = "S10.2", name = "Copy to Publication Directory", target = "prot-summary-copy-publication",
            dom = \(st) nzchar(k233Text(st, copy_id)) && !grepl("^Error", k233Text(st, copy_id)) &&
                grepl("Files copied", k233Text(st, summary_id), fixed = TRUE),
            outputs = c(copy_id, summary_id)),
        list(step = "S10.3", name = "Generate Report", target = "prot-summary-generate-report",
            sels = "[data-testid='prot-summary-download-report']",
            dom = \(st) isTRUE(st$vis[[1L]]) || k233Noted(st, "Report generated successfully", "message"),
            outputs = c(summary_id, sm("template_status"))),
        list(step = "S10.4", name = "Download Report", custom = function(S) k233DownloadReport(S, R)),
        list(step = "S10.5", name = "Export Session State", target = "prot-summary-export-state",
            dom = \(st) k233Noted(st, "Session state exported to:", "message"))
    )
    res <- list(outcome = "complete")
    for (a in actions) {
        a$by <- "testid"
        a$tabs <- tabs
        a$texts <- unique(c(summary_id, copy_id))
        res <- k233Action(S, a)
        if (!identical(res$outcome, "complete")) return(res)
    }
    res
}

k233DownloadReport <- function(S, R) {
    p <- k233Press(S, "S10.4", "prot-summary-download-report", "testid")
    if (!isTRUE(p$ok)) return(list(outcome = "refused", message = paste("not reachable as a user:", p$reason)))
    path <- tryCatch(S$G$helpers$e2e_get_download(S$driver, paste0(K233_NS, "session_summary-download_report"),
        artifact_dir = S$G$work, filename = paste0(S$label, "-report-download")), error = \(e) conditionMessage(e))
    if (!file.exists(path) || file.size(path) == 0) {
        return(list(outcome = "refused", message = paste("the download gave no file:", path)))
    }
    magic <- readBin(path, "raw", 5L)
    type <- if (identical(magic[1:2], charToRaw("PK"))) "docx" else if (identical(magic[1:4], charToRaw("%PDF"))) "pdf" else
        "html"
    S$report <- list(path = path, type = type, bytes = file.size(path), sha256 = k233Sha(path))
    text <- switch(type,
        docx = {
            dir <- file.path(S$G$work, paste0(S$label, "-report-docx"))
            utils::unzip(path, files = "word/document.xml", exdir = dir)
            xml <- paste(readLines(file.path(dir, "word/document.xml"), warn = FALSE), collapse = "")
            gsub("<[^>]+>", "", gsub("</w:p>", "\n", xml))
        },
        html = gsub("<[^>]+>", " ", paste(readLines(path, warn = FALSE), collapse = "\n")),
        "")
    if (identical(type, "html")) file.copy(path, file.path(S$G$out, paste0(S$label, "-report.html")), overwrite = TRUE)
    writeLines(text, file.path(S$G$out, paste0(S$label, "-report-text.txt")))
    S$report$text <- text
    k233NumberN4(S, R)
    list(outcome = "complete", message = "", dom_t = NA, dig_t = NA, agreement = "dom_only",
        note = sprintf("report %s, %d bytes, sha256 %s", type, S$report$bytes, S$report$sha256))
}

# --- the interface numbers against the saved files (5.4) --------------------------------------------------------

k233Number <- function(S, id, contrast, shown, saved, where, planted = FALSE) {
    result <- if (is.na(shown) || is.na(saved)) "not_computed" else if (identical(as.character(shown),
        as.character(saved))) "equal" else "differ"
    S$G$numbers[[length(S$G$numbers) + 1L]] <- list(session = S$label, id = id, contrast = contrast, shown = shown,
        saved = saved, rule = "equal", result = result, where = where, planted = planted)
    if (identical(result, "differ")) {
        k233Finding(S, switch(id, N1 = "S7", N2 = "S8.3", N3 = "S10.1", N4 = "S10.4", id), 4L,
            sprintf("%s %s: shown %s, saved %s (%s)", id, contrast, shown, saved, where), planted = planted)
    }
    result
}

k233NumberN1 <- function(S, R) {
    file <- file.path(k233ProjectPaths(R)$source, "filtered_session_summary.txt")
    lines <- k233ReadLog(file)
    field <- function(label) {
        hit <- grep(paste0("^- ", label, ": "), lines, value = TRUE)
        if (length(hit)) suppressWarnings(as.integer(sub(paste0("^- ", label, ": "), "", hit[[1L]]))) else NA_integer_
    }
    if (file.exists(file)) file.copy(file, file.path(S$G$out, paste0(S$label, "-filtered_session_summary.txt")),
        overwrite = TRUE)
    shown <- S$n1_shown %||% list(proteins = NA_integer_, samples = NA_integer_)
    where <- "normalization-correlation_filter_summary vs scripts/proteomics/filtered_session_summary.txt"
    k233Number(S, "N1", "proteins", shown$proteins, field("Proteins"), where)
    k233Number(S, "N1", "samples", shown$samples, field("Samples"), where)
}

k233DaCount <- function(path, threshold, lfc) {
    table <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    qcol <- intersect(c("fdr_qvalue", "q.mod", "adj.P.Val", "qvalue"), names(table))
    if (!length(qcol)) return(list(count = NA_integer_, rule = "no q column"))
    q <- suppressWarnings(as.numeric(table[[qcol[[1L]]]]))
    keep <- !is.na(q) & q < threshold
    if (isTRUE(lfc > 0) && "log2FC" %in% names(table)) keep <- keep & abs(as.numeric(table$log2FC)) >= lfc
    list(count = sum(keep), rule = sprintf("%s < %s%s", qcol[[1L]], threshold, if (isTRUE(lfc > 0))
        sprintf(" and |log2FC| >= %s", lfc) else ""), column = qcol[[1L]])
}

k233DaFile <- function(R, contrast) {
    dir <- k233ProjectPaths(R)$da
    files <- list.files(dir, "_long_annot\\.tsv$", full.names = TRUE)
    exact <- file.path(dir, paste0("da_proteins_", gsub("[^A-Za-z0-9_-]", "_", contrast), "_long_annot.tsv"))
    if (file.exists(exact)) return(exact)
    hit <- files[k233Norm(sub("^da_proteins_(.*)_long_annot\\.tsv$", "\\1", basename(files))) == k233Norm(contrast)]
    if (length(hit)) hit[[1L]] else NA_character_
}

k233NumberN2 <- function(S, R) {
    da <- function(x) paste0(K233_NS, "differential_expression-", x)
    threshold <- suppressWarnings(as.numeric(k233InputValue(S, da("da_q_val_thresh"))))
    lfc <- suppressWarnings(as.numeric(k233InputValue(S, da("treat_lfc_cutoff"))))
    options <- as.character(unlist(k233Js(S, sprintf(
        "Array.prototype.map.call((document.getElementById(%s)||{options:[]}).options,function(o){return o.value;})",
        k233Json(da("table_contrast"))))))
    S$n2 <- list()
    for (contrast in options) {
        k233Set(S, "S8.3", da("table_contrast"), contrast)
        started <- Sys.time()
        text <- ""
        repeat {
            text <- as.character(k233Js(S, sprintf("(document.getElementById(%s)||{innerText:''}).innerText",
                k233Json(da("da_summary_stats")))) %||% "")
            if (grepl("Significant", text, fixed = TRUE) || k233Seconds(started) > 60) break
            Sys.sleep(1)
        }
        Sys.sleep(1)
        text <- as.character(k233Js(S, sprintf("(document.getElementById(%s)||{innerText:''}).innerText",
            k233Json(da("da_summary_stats")))) %||% text)
        shown <- suppressWarnings(as.integer(sub("(?s).*Significant \\(q < [0-9.]+\\): ([0-9]+).*", "\\1", text, perl = TRUE)))
        file <- k233DaFile(R, contrast)
        saved <- if (is.na(file)) list(count = NA_integer_, rule = "no saved DA table") else k233DaCount(file, threshold, lfc)
        if (!is.na(file)) file.copy(file, file.path(S$G$out, paste0(S$label, "-", basename(file))), overwrite = TRUE)
        k233Number(S, "N2", contrast, shown, saved$count, sprintf("da_summary_stats vs %s (%s)",
            if (is.na(file)) "no file" else basename(file), saved$rule))
        S$n2[[contrast]] <- list(shown = shown, file = file, threshold = threshold, lfc = lfc, saved = saved$count)
    }
    invisible(S$n2)
}

k233NumberN3 <- function(S, R) {
    sm <- function(x) paste0(K233_NS, "session_summary-", x)
    file <- file.path(k233ProjectPaths(R)$source, "study_parameters.txt")
    lines <- k233ReadLog(file)
    field <- function(label) {
        hit <- grep(paste0("^", label, ": "), lines, value = TRUE)
        if (length(hit)) sub(paste0("^", label, ": "), "", hit[[1L]]) else NA_character_
    }
    if (file.exists(file)) file.copy(file, file.path(S$G$out, paste0(S$label, "-study_parameters.txt")), overwrite = TRUE)
    where <- "session_summary inputs vs scripts/proteomics/study_parameters.txt"
    k233Number(S, "N3", "label", k233Clean(k233InputValue(S, sm("experiment_label"))), field("Workflow Name"), where)
    k233Number(S, "N3", "description", k233Clean(k233InputValue(S, sm("description"))), field("Description"), where)
}

k233NumberN4 <- function(S, R) {
    text <- gsub("\\s+", " ", S$report$text %||% "")
    pattern <- "For the comparison (.+?), there were ([0-9]+) upregulated and ([0-9]+) downregulated proteins\\."
    hits <- regmatches(text, gregexpr(pattern, text, perl = TRUE))[[1L]]
    none <- grepl("No significant differentially expressed proteins were found", text, fixed = TRUE)
    named <- k233Norm(sub(pattern, "\\1", hits, perl = TRUE))
    same <- function(a, b) a == b || grepl(a, b, fixed = TRUE) || grepl(b, a, fixed = TRUE)
    contrasts <- names(S$n2 %||% list())
    # The report names a comparison only when it has significant proteins; an unnamed comparison counts 0 only when
    # the report's names are in the tab's form (at least one name matches a contrast), else it is not compared.
    compatible <- any(vapply(named, \(n) any(vapply(k233Norm(contrasts), \(c) same(n, c), logical(1))), logical(1)))
    for (contrast in contrasts) {
        n2 <- S$n2[[contrast]]
        shown <- NA_integer_
        for (i in seq_along(hits)) {
            if (same(named[[i]], k233Norm(contrast))) {
                shown <- as.integer(sub(pattern, "\\2", hits[[i]], perl = TRUE)) +
                    as.integer(sub(pattern, "\\3", hits[[i]], perl = TRUE))
                break
            }
        }
        if (is.na(shown) && (none || compatible)) shown <- 0L
        k233Number(S, "N4", contrast, shown, n2$shown, sprintf("report sentence (%s) vs N2 shown count%s",
            S$report$type %||% "no report", if (length(hits) && !compatible) "; report names not in the tab's form" else ""))
    }
}

# ---------------------------------------------------------------------------------------------------------------
# The runs.
# ---------------------------------------------------------------------------------------------------------------
k233FixtureRun <- function(G, project_dir) {
    dir <- file.path(G$pkg, "tests/testdata/e2e/prot_dia")
    design <- utils::read.delim(file.path(dir, "design.tsv"), stringsAsFactors = FALSE)
    groups <- unique(design$group)
    list(dataset = "prot_dia", report = file.path(dir, "report.tsv"), fasta = file.path(dir, "proteins.fasta"),
        groups = groups, members = lapply(stats::setNames(groups, groups), \(g) design$sample[design$group == g]),
        contrasts = list(c("KO", "WT")), label = "k233_lifecycle",
        description = "KIT-233 lifecycle on the prot_dia fixture", permissive = TRUE, project_dir = project_dir)
}

k233PxdRun <- function(G, project_dir) {
    ins <- G$phase$inputs
    groups <- sub("-", "_", K233_MIXES, fixed = TRUE)
    members <- lapply(stats::setNames(K233_MIXES, groups), \(m) ins$runs[ins$mix == m])
    list(dataset = "PXD055981 s500", report = ins$report, fasta = ins$fasta, groups = groups, members = members,
        contrasts = list(c("H1E_E", "H1E_C"), c("H1E_C", "H1E_B"), c("H1E_D", "H1E_B"), c("H1E_E", "H1E_B"),
            c("H1E_D", "H1E_C"), c("H1E_E", "H1E_D")),
        label = "k233_s500", description = "KIT-233 discovery journey on PXD055981, 500 SPD 15 ng, mixes H1E-B to E",
        permissive = FALSE, project_dir = project_dir)
}

k233Journey <- function(G) {
    S <- k233Launch(G, "J1")
    if (!isTRUE(S$answered) || is.null(S$driver)) {
        writeLines(sprintf("instrument: the application did not answer or could not be attached (%s)", S$attach_error),
            file.path(G$out, "verdict.txt"))
        return(invisible(NULL))
    }
    k233RunRoute(S, G$route, k233PxdRun(G, S$project_dir))
    k233StopSession(S)
    writeLines("journey recorded (a result run: findings in findings.tsv)", file.path(G$out, "verdict.txt"))
}

k233L0 <- function(G, label = "L0", test_mode = TRUE) {
    crit <- function(item, pass, evidence = "") k233Criterion(G, "L-0", item, pass, evidence)
    chrome <- tryCatch(chromote::find_chrome(), error = \(e) "")
    crit("Chrome found", nzchar(chrome) && file.exists(chrome), chrome)
    crit("shinytest2 loads", requireNamespace("shinytest2", quietly = TRUE),
        tryCatch(as.character(utils::packageVersion("shinytest2")), error = \(e) "not installed"))
    S <- k233Launch(G, label, test_mode = test_mode)
    log <- k233ReadLog(S$log)
    crit("the entry.R launch passes its trailing arguments, starts and answers on the port",
        isTRUE(S$answered) && !isTRUE(S$usage_error), sprintf("answered=%s after %s s; usage error=%s; %s", S$answered,
            S$answer_s, S$usage_error, S$command))
    crit("the MSR-LAUNCHER line names process", grepl('"launcher":"process"', S$launcher, fixed = TRUE), S$launcher)
    crit("legacy RUV route: no MSR-TEST-FEATURES line and MULTISCHOLAR_TEST_FEATURES unset (the panel's absence in the page is checked at S6 of L-1)",
        !any(grepl("^MSR-TEST-FEATURES", log)) && !nzchar(Sys.getenv("MULTISCHOLAR_TEST_FEATURES")),
        "entry.R sets MultiScholaR.ruv_design_panel only under MULTISCHOLAR_TEST_FEATURES=1 (packaging/electron/entry.R:20-22)")
    crit("AppDriver attaches by URL", !is.null(S$driver), S$attach_error)
    if (is.null(S$driver)) return(invisible(S))
    st <- k233State(S, 0)
    if (test_mode) {
        crit("test_state_digest is read; the controller sees test mode", !is.null(st$dg),
            if (is.null(st$dg)) "no digest in the page" else sprintf("digest keys: %s", paste(names(st$dg), collapse = ",")))
    } else {
        crit("user mode: no digest in the page", is.null(st$dg), "the digest exists only in test mode (R/utils_test_mode.R:71)")
    }
    crit("no worker code reads test mode", NA, paste(
        "not executable by the driver; code read: worker publication forces the option off",
        "(R/utils_workflow_execution.R:3149-3150, :3214-3215), but workers inherit the controller's environment",
        "(processx env 'current', :508) and is_test_mode() also reads MULTISCHOLAR_TEST_MODE (R/utils_test_mode.R:22-27)"))
    crit(sprintf("network: %s", G$phase$launch$network), G$phase$launch$network %in% c("namespace", "proxy", "firewall"),
        k233Clean(vapply(G$env, \(r) paste(r$key, r$value, sep = "="), character(1))[2:4]))
    crit("the inner phase's effective user equals the outer phase's", identical(Sys.info()[["effective_user"]],
        G$phase$user_outer), sprintf("inner %s, outer %s", Sys.info()[["effective_user"]], G$phase$user_outer))
    if (!test_mode) {
        return(invisible(S))
    }
    R <- k233FixtureRun(G, S$project_dir)
    R$label <- "k233_upload_check"
    home <- k233StepHome(S, R, list(id = "S1", name = "home and project (for the upload check)"), G$route)
    big <- file.path(G$work, "k233-20mb-report.tsv")
    lines <- readLines(R$report)
    body <- lines[-1L]
    times <- ceiling(20 * 1024^2 / sum(nchar(body, type = "bytes") + 1)) + 1L
    writeLines(c(lines[[1L]], rep(body, times)), big)
    up <- if (identical(home$outcome, "complete")) {
        k233ActivateTab(S, k233Workflow("setup")[[1L]], "setup")
        k233Upload(S, "L0", paste0(K233_NS, "setup_import-search_results_standard"), big, timeout_s = 240)
    } else {
        list(ok = FALSE, text = paste("S1 was", home$outcome))
    }
    crit("a 20 MB DIA-NN table uploads with no size error (the raised limit)",
        isTRUE(up$ok) && !grepl("Maximum upload size", up$text, ignore.case = TRUE),
        sprintf("%.1f MB; progress text: %s; %s s", file.size(big) / 1024^2, up$text, up$seconds %||% ""))
    Sys.sleep(5)
    k233Observation(G, label, "format detection after the 20 MB upload", k233Clean(k233Js(S, sprintf(
        "(document.getElementById(%s)||{innerText:''}).innerText", k233Json(paste0(K233_NS, "setup_import-format_detection"))))))
    k233StopSession(S)
    invisible(S)
}

k233L1Criteria <- function(G, label = "L1", criterion = "L-1") {
    crit <- function(item, pass, evidence = "") k233Criterion(G, criterion, item, pass, evidence)
    rows <- Filter(\(r) identical(r$session, label), G$steps)
    S <- G$sessions[[label]]
    outcome <- vapply(rows, \(r) r$outcome, character(1))
    step <- vapply(rows, \(r) r$step, character(1))
    crit("every step's row has an outcome and seconds", length(rows) > 0 && all(nzchar(outcome)) &&
        all(vapply(rows, \(r) !is.null(r$seconds) && !is.na(r$seconds), logical(1))), paste(step, outcome, collapse = "; "))
    crit("reached S10's last action", any(step == "S10.5" & outcome == "complete"),
        paste(utils::tail(paste(step, outcome), 3L), collapse = "; "))
    s3 <- Filter(\(r) identical(r$step, "S3"), rows)
    if (S$test_mode) {
        crit("S3 complete with KIT-175's annotation warning in its notifications", length(s3) == 1L &&
            identical(s3[[1L]]$outcome, "complete") && grepl("UniProt|annotation", s3[[1L]]$warnings, ignore.case = TRUE),
            if (length(s3)) paste(s3[[1L]]$outcome, s3[[1L]]$warnings) else "no S3 row")
    }
    crit("the organism modal answered and recorded at S2", isTRUE(S$organism_answered), k233Clean(S$organism_modal))
    complete <- Filter(\(r) identical(r$outcome, "complete"), rows)
    listed <- unlist(lapply(complete, \(r) sub(" .*$", "", strsplit(k233Clean(r$files), " \\| ")[[1L]])))
    present <- file.exists(file.path(G$out, listed)) & file.size(file.path(G$out, listed)) > 0
    kinds <- vapply(complete, \(r) all(vapply(c("\\.png ", "notifications\\.json ", "console\\.json ", "outputs\\.json "),
        \(p) grepl(p, paste0(k233Clean(r$files), " ")), logical(1))), logical(1))
    crit("every screenshot and JSON exists and is non-empty", length(listed) > 0 && all(present) && all(kinds),
        sprintf("%d files, %d present; rows with all four kinds: %d of %d", length(listed), sum(present), sum(kinds),
            length(kinds)))
    crit("the R log is present", file.exists(S$log) && file.size(S$log) > 0, sprintf("%s, %s bytes", basename(S$log),
        if (file.exists(S$log)) file.size(S$log) else 0))
    numbers <- Filter(\(n) identical(n$session, label), G$numbers)
    ids <- vapply(numbers, \(n) n$id, character(1))
    computed <- vapply(numbers, \(n) !identical(n$result, "not_computed"), logical(1))
    crit("N1 to N4 computed", all(c("N1", "N2", "N3", "N4") %in% ids[computed]) && all(computed),
        paste(ids, vapply(numbers, \(n) n$result, character(1)), collapse = "; "))
    pressed <- Filter(\(r) !identical(r$outcome, "not_reached") && !startsWith(r$outcome, "not_reached"), rows)
    checks <- vapply(pressed, \(r) k233Clean(r$tab_checks), character(1))
    crit("the active-tab check recorded for every press", all(nzchar(checks)) && !any(grepl("NOT ", checks, fixed = TRUE)),
        paste(utils::head(checks[grepl("NOT ", checks, fixed = TRUE)], 3L), collapse = " || "))
    progress <- vapply(rows, \(r) k233Clean(r$progress_modals), character(1))
    terminal <- vapply(rows, \(r) identical(r$outcome, "refused") && startsWith(k233Clean(r$message), "modal '"), logical(1))
    crit("every progress modal seen recorded as such, none classified terminal", !any(terminal),
        paste(progress[nzchar(progress)], collapse = " || "))
    invisible(NULL)
}

k233L4a <- function(S, R) {
    da <- paste0(K233_NS, "differential_expression-da_results_table")
    emptied <- k233Js(S, sprintf(paste0("(function(){var e=document.getElementById(%s);if(!e){return false;}",
        "var b=e.querySelectorAll('tbody');b.forEach(function(x){x.innerHTML='';});return b.length>0;})()"), k233Json(da)))
    outputs <- k233CheckOutputs(S, c(paste0(K233_NS, "differential_expression-da_summary_stats"), da))
    empty <- vapply(Filter(\(o) !isTRUE(o$nonempty), outputs), \(o) o$id, character(1))
    for (o in empty) k233Finding(S, "S8", 3L, paste("named output empty:", o, "(L-4 (a), planted)"), planted = TRUE)
    k233Criterion(S$G, "L-4", "(a) the emptied DA results table gives a kind-3 finding naming it", da %in% empty,
        sprintf("tbody emptied: %s; empty outputs at the second check: %s", k233Clean(emptied), paste(empty, collapse = ", ")))
}

k233L4b <- function(G, S) {
    n2 <- Filter(\(x) !is.na(x$file), S$n2 %||% list())
    if (!length(n2)) return(k233Criterion(G, "L-4", "(b) a q moved across 0.05 gives a kind-4 finding", FALSE,
        "no saved DA table from L-1"))
    first <- n2[[1L]]
    table <- utils::read.delim(first$file, check.names = FALSE, stringsAsFactors = FALSE)
    column <- k233DaCount(first$file, first$threshold, first$lfc)$column
    q <- suppressWarnings(as.numeric(table[[column]]))
    below <- which(!is.na(q) & q < first$threshold)
    row <- if (length(below)) below[[1L]] else which.max(q)
    table[[column]][row] <- if (length(below)) first$threshold + 0.01 else first$threshold / 2
    copy <- file.path(G$work, paste0("planted-", basename(first$file)))
    utils::write.table(table, copy, sep = "\t", quote = FALSE, row.names = FALSE)
    planted <- k233DaCount(copy, first$threshold, first$lfc)$count
    result <- k233Number(S, "N2", paste(names(n2)[[1L]], "(planted copy)"), first$shown, planted,
        sprintf("da_summary_stats vs a copy of %s with row %d's q moved across %s", basename(first$file), row,
            first$threshold), planted = TRUE)
    k233Criterion(G, "L-4", "(b) a q moved across 0.05 gives a kind-4 finding", identical(result, "differ"),
        sprintf("shown %s, planted copy %s", first$shown, planted))
}

k233L5 <- function(G, label = "L1") {
    rows <- Filter(\(r) identical(r$session, label) && identical(r$outcome, "complete"), G$steps)
    agreement <- vapply(rows, \(r) k233Clean(r$agreement), character(1))
    bad <- vapply(rows, \(r) r$step, character(1))[!agreement %in% c("agree", "dom_only")]
    k233Criterion(G, "L-5", "the digest and the DOM signals agree on completion at every step of L-1", !length(bad),
        if (length(bad)) paste("disagree at", paste(bad, collapse = ", ")) else
            sprintf("%d agree, %d DOM only", sum(agreement == "agree"), sum(agreement == "dom_only")))
}

k233KillAppHook <- function(S, el) {
    if (isTRUE(S$killed)) return(invisible(NULL))
    S$killed <- TRUE
    S$kill_at <- el
    S$proc$kill()
    invisible(NULL)
}

k233KillWorkerHook <- function(S, el) {
    if (isTRUE(S$worker_killed) || el > 120) return(invisible(NULL))
    kids <- tryCatch(ps::ps_children(ps::ps_handle(S$proc$get_pid()), recursive = TRUE), error = \(e) list())
    for (k in kids) {
        command <- tryCatch(paste(ps::ps_cmdline(k), collapse = " "), error = \(e) "")
        if (grepl("workflowWorkerMain", command, fixed = TRUE)) {
            tryCatch(ps::ps_kill(k), error = \(e) NULL)
            S$worker_killed <- TRUE
            S$worker_kill_at <- round(el, 1)
            break
        }
    }
    invisible(NULL)
}

k233Lifecycle <- function(G) {
    S0 <- k233L0(G)
    if (is.null(S0$driver)) {
        writeLines("instrument not proven: L-0 (launch or attach) failed", file.path(G$out, "verdict.txt"))
        return(invisible(NULL))
    }
    S1 <- k233Launch(G, "L1")
    R <- k233FixtureRun(G, S1$project_dir)
    R$after <- list(S8 = k233L4a)
    k233RunRoute(S1, G$route, R)
    k233StopSession(S1)
    guard <- function(criterion, code) tryCatch(code, error = \(e) k233Criterion(G, criterion,
        "driver error while judging", FALSE, conditionMessage(e)))
    guard("L-1", k233L1Criteria(G, "L1"))
    cover <- Filter(\(c) identical(c$session, "L1"), G$coverage)
    tall <- Filter(\(c) isTRUE(suppressWarnings(as.numeric(c$pane_height)) > 1100), cover)
    k233Criterion(G, "L-0", "for a tall tab, each named output's box lies inside a saved image (pane or its own)",
        length(tall) > 0 && all(vapply(tall, \(c) !identical(c$covered, "not covered"), logical(1))),
        sprintf("%d outputs on panes taller than the viewport; not covered: %s", length(tall),
            paste(vapply(Filter(\(c) identical(c$covered, "not covered"), cover), \(c) c$output, character(1)),
                collapse = ", ")))
    k233Criterion(G, "L-0", "the RUV design panel is absent at S6 (legacy route)", !isTRUE(G$ruv_panel),
        sprintf("ruv_design_summary present: %s", k233Clean(G$ruv_panel)))
    guard("L-4", k233L4b(G, S1))
    guard("L-5", k233L5(G, "L1"))
    if (!any(vapply(G$criteria, \(c) identical(c$criterion, "L-4") && startsWith(c$item, "(a)"), logical(1)))) {
        k233Criterion(G, "L-4", "(a) the emptied DA results table gives a kind-3 finding naming it", FALSE,
            "S8 of L-1 did not complete, so the planted check did not run")
    }
    S2 <- k233Launch(G, "L2")
    missing <- file.path(G$work, "no-such-folder")
    R2 <- R
    R2$project_dir <- missing
    k233RunRoute(S2, G$route, R2)
    k233StopSession(S2)
    rows <- Filter(\(r) identical(r$session, "L2"), G$steps)
    s1 <- Filter(\(r) identical(r$step, "S1"), rows)
    later <- Filter(\(r) !identical(r$step, "S1"), rows)
    note_t <- NA_real_
    notes_file <- file.path(G$out, "L2-01-S1-notifications.json")
    if (file.exists(notes_file)) {
        during <- jsonlite::read_json(notes_file)$during_step
        error <- Filter(\(n) identical(n$type, "error"), during)
        if (length(error)) note_t <- error[[1L]]$t
    }
    press <- Filter(\(a) identical(a$session, "L2") && identical(a$input, "btn-confirm-setup"), G$actions)
    k233Criterion(G, "L-2", "S1 refused with 'Base directory does not exist: <path>'; S2 to S10 not_reached",
        length(s1) == 1L && identical(s1[[1L]]$outcome, "refused") &&
            grepl(paste("Base directory does not exist:", missing), s1[[1L]]$message, fixed = TRUE) &&
            length(later) == 9L && all(vapply(later, \(r) identical(r$outcome, "not_reached"), logical(1))),
        sprintf("S1 %s: %s; later rows: %s; error notification at browser time %s", if (length(s1)) s1[[1L]]$outcome else "-",
            if (length(s1)) s1[[1L]]$message else "-", paste(vapply(later, \(r) r$outcome, character(1)), collapse = ","),
            note_t))
    S3 <- k233Launch(G, "L3")
    R3 <- R
    R3$project_dir <- S3$project_dir
    k233RunRoute(S3, G$route, R3, upto = "S2", hooks = list(S2 = k233KillAppHook))
    k233StopSession(S3)
    s2 <- Filter(\(r) identical(r$session, "L3") && identical(r$step, "S2"), G$steps)
    k233Criterion(G, "L-3", "the controller killed while S2 waits: S2 recorded crashed (not timeout)",
        length(s2) == 1L && identical(s2[[1L]]$outcome, "crashed"),
        sprintf("S2 %s: %s; killed %s s after the press", if (length(s2)) s2[[1L]]$outcome else "-",
            if (length(s2)) s2[[1L]]$message else "-", k233Clean(S3$kill_at)))
    k233LifecycleVerdict(G)
    if (k233Remaining(G) > 300) {
        So <- k233Launch(G, "L2o")
        Ro <- R
        Ro$project_dir <- So$project_dir
        Ro$report <- file.path(G$pkg, "tests/testdata/e2e/prot_dia/design.tsv")
        k233RunRoute(So, G$route, Ro, upto = "S2")
        k233StopSession(So)
        k233Observation(G, "L2o", "design.tsv given as the search results (L-2 observation)", k233Clean(vapply(
            Filter(\(r) identical(r$session, "L2o"), G$steps), \(r) paste(r$step, r$outcome, r$message), character(1))))
    }
    qc <- Filter(\(r) identical(r$session, "L1") && grepl("^S[45]\\.[0-9]+$", r$step) && identical(r$outcome, "complete"),
        G$steps)
    if (length(qc) && k233Remaining(G) > 420) {
        longest <- qc[[which.max(vapply(qc, \(r) as.numeric(r$seconds), numeric(1)))]]
        Sw <- k233Launch(G, "L3w")
        Rw <- R
        Rw$project_dir <- Sw$project_dir
        hooks <- stats::setNames(list(k233KillWorkerHook), longest$step)
        k233RunRoute(Sw, G$route, Rw, upto = sub("\\..*$", "", longest$step), hooks = hooks)
        k233StopSession(Sw)
        shown <- Filter(\(r) identical(r$session, "L3w") && identical(r$step, longest$step), G$steps)
        k233Observation(G, "L3w", sprintf("worker killed during %s (%s, the longest QC step of L-1, %s s)", longest$step,
            longest$name, longest$seconds), sprintf("worker killed: %s at %s s; the step: %s", isTRUE(Sw$worker_killed),
            k233Clean(Sw$worker_kill_at), if (length(shown)) paste(shown[[1L]]$outcome, shown[[1L]]$message) else "-"))
    }
    k233LifecycleVerdict(G)
}

k233LifecycleVerdict <- function(G) {
    results <- vapply(G$criteria, \(c) c$result, character(1))
    names <- vapply(G$criteria, \(c) c$criterion, character(1))
    blocking <- names %in% c("L-0", "L-1", "L-2", "L-3", "L-4")
    proven <- !any(results[blocking] == "fail")
    l5 <- !any(results[names == "L-5"] == "fail")
    writeLines(c(sprintf("instrument %s (L-0 to L-4: %d pass, %d fail, %d recorded)", if (proven) "proven" else "not proven",
        sum(results[blocking] == "pass"), sum(results[blocking] == "fail"), sum(results[blocking] == "recorded")),
        sprintf("L-5 %s", if (l5) "passed" else "not passed: J may go; no user-mode run may rely on the DOM signals")),
        file.path(G$out, "verdict.txt"))
    k233WriteAll(G)
}

#' L-u (observation only): L-0 and L-1 in user mode, test mode off, the real file chooser.
k233UserLife <- function(G) {
    S0 <- k233L0(G, "U0", test_mode = FALSE)
    k233StopSession(S0)
    S1 <- k233Launch(G, "U1", test_mode = FALSE)
    if (is.null(S1$driver)) return(invisible(NULL))
    R <- k233FixtureRun(G, S1$project_dir)
    home <- G$phase$inputs$home_dir
    R$report <- file.path(home, "report.tsv")
    R$fasta <- file.path(home, "proteins.fasta")
    R$label <- "k233_user"
    k233RunRoute(S1, G$route, R)
    k233StopSession(S1)
    k233L1Criteria(G, "U1", "L-u")
    writeLines("L-u recorded (observation only)", file.path(G$out, "verdict.txt"))
    k233WriteAll(G)
}

if (sys.nframe() == 0L) k233Main(commandArgs(trailingOnly = TRUE))
