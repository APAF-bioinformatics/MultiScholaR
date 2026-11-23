# KIT-180 script S2: the control and the freeze of the GUI-entry references prot_tmt, prot_lfq and prot_lfq_fragpipe
# (ticket section 6.2; KIT-159 r1.5 gate G1(a), chain link (i) and conditions (ii) and (iii)).
#
# Usage, from the package root of the fix commit on a run-tests runner (full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-180/gui-entry-control.R <base commit> [<out dir>]
#   <out dir>  default docs/qa/workflow-execution/kit/kit-180/control (made absolute)
#
# (a) Base side: a git worktree of <base commit>; a fresh Rscript process with its working directory at that tree's root
#     loads it (pkgload::load_all, helpers sourced into the global environment) and runs guiEntryRunLane() for the three
#     lanes with MULTISCHOLAR_GUI_ENTRY_RECORD_DIR=<out>/base. No comparison is made there. Each lane runs inside
#     testthat::test_that(), as the head side and the test files do: R4 (run 45129403566) showed that outside it the
#     store-entry names of a capture are sorted in the runner's locale, not in test_that()'s C collation, so the base
#     record differed from the frozen files and the head in key order only. R4 is a control result that stopped; the
#     run with this change (R4b) is the further control run the owner ratified (VAL-63).
# (b) Head side: a fresh Rscript process at this tree's root runs guiEntryExpectLane() for each lane inside
#     testthat::test_that() with MULTISCHOLAR_GUI_ENTRY_RECORD_DIR=<out>/head, MULTISCHOLAR_GUI_ENTRY_FREEZE=1 and
#     MULTISCHOLAR_GUI_ENTRY_ORACLE_DIR=<out>/oracles: one run per lane writes its record and its frozen file into out,
#     never into the checkout.
# (c) Comparison, here, from the base tree's frozen files, <out>/base, <out>/head and <out>/oracles, one row per differing
#     key (packaging/verify/run-lanes.R compareLane's form), each with its class:
#       frozen_vs_base   step 1: G1(a)'s classes by their bounds (integration_artifacts: `artifacts` keys
#                        integration/proteomics/... present on one side only, and `<store entries>`); anything else is
#                        `outside`. Per lane, the steps of each class against kit-180/g1a-reference-protein-lanes.tsv.
#       base_vs_head     step 2: the changed step is `design` on each lane; from it on, `tab_status` quality_control and
#                        normalization and the stepper output proteomics_workflow-workflow_progress at the five steps
#                        design, qc_protein_cleanup, qc_protein_intensity_filter_refusal, qc_protein_intensity_filter,
#                        qc_protein_duplicate_removal, with head quality_control `pending` and normalization `disabled`,
#                        are `qc_routing (KIT-180)`; anything else is `outside` (before `design`, every capture equal).
#       head_vs_oracle   the head record against the file frozen from it (expected: none).
#     Written: <out>/gui-entry-control.tsv (comparison, lane, step, field, key, kind, a, b, class, changed_step,
#     routing_step), <out>/gui-entry-steps.tsv (lane, step, frozen, base and head step digests, cause, routing_step) and
#     <out>/verdict.tsv (check, value). Only <out>/oracles/*.json are committed, to tests/testdata/gui-entry/, and only
#     when scripts.tsv shows this script `ok` and the verdict holds.
# Exit status non-zero when a head-side test_that() block returns FALSE or a lane errors on either side.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L, grepl("^[0-9a-f]{40}$", args[[1L]]))
base_commit <- args[[1L]]
out <- if (length(args) >= 2L) args[[2L]] else "docs/qa/workflow-execution/kit/kit-180/control"
lanes <- c("prot_tmt", "prot_lfq", "prot_lfq_fragpipe")
routing_steps <- c("design", "qc_protein_cleanup", "qc_protein_intensity_filter_refusal", "qc_protein_intensity_filter",
    "qc_protein_duplicate_removal")
for (sub in c("base", "head", "oracles")) dir.create(file.path(out, sub), recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
head_root <- normalizePath(".")
head_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
cat("KIT-180 S2: head", head_commit, "base", base_commit, "out", out, "\n")
rscript <- file.path(R.home("bin"), "Rscript")

child <- file.path(tempdir(), "kit180-s2-child.R")
writeLines(c(
    "args <- commandArgs(trailingOnly = TRUE)",
    "setwd(args[[1L]])",
    "side <- args[[2L]]",
    "lanes <- strsplit(args[[3L]], ',', fixed = TRUE)[[1L]]",
    "suppressPackageStartupMessages(pkgload::load_all('.', quiet = TRUE, export_all = TRUE, helpers = FALSE))",
    "suppressMessages(testthat::source_test_helpers('tests/testthat', env = globalenv()))",
    "failed <- FALSE",
    "for (lane in lanes) {",
    "    ok <- tryCatch(if (identical(side, 'head')) {",
    "        testthat::test_that(paste('KIT-180 control', lane), guiEntryExpectLane(lane))",
    "    } else {",
    "        testthat::test_that(paste('KIT-180 control base', lane), { guiEntryRunLane(lane); testthat::succeed() })",
    "    }, error = function(e) {",
    "        cat('lane', lane, 'error', class(e)[[1L]], conditionMessage(e), '\\n')",
    "        FALSE",
    "    })",
    "    cat('KIT-180 S2', side, lane, if (isTRUE(ok)) 'ok' else 'FAILED', '\\n')",
    "    if (!isTRUE(ok)) failed <- TRUE",
    "}",
    "quit(status = if (failed) 1L else 0L)"), child)
runSide <- function(tree, side, env) {
    console <- suppressWarnings(system2(rscript, c(shQuote(child), shQuote(tree), side, paste(lanes, collapse = ",")),
        stdout = TRUE, stderr = TRUE, env = env))
    cat(grep("^KIT-180 S2|^lane ", console, value = TRUE), sep = "\n")
    cat(utils::tail(console, 10L), sep = "\n")
    attr(console, "status") %||% 0L
}

# (a) Base side.
base_tree <- file.path(tempdir(), "kit180-s2-base")
invisible(system2("git", c("worktree", "add", "--detach", shQuote(base_tree), base_commit), stdout = TRUE, stderr = TRUE))
stopifnot(dir.exists(file.path(base_tree, "tests", "testthat")))
base_status <- runSide(base_tree, "base", paste0("MULTISCHOLAR_GUI_ENTRY_RECORD_DIR=", shQuote(file.path(out, "base"))))

# (b) Head side.
head_status <- runSide(head_root, "head", c(
    paste0("MULTISCHOLAR_GUI_ENTRY_RECORD_DIR=", shQuote(file.path(out, "head"))),
    "MULTISCHOLAR_GUI_ENTRY_FREEZE=1",
    paste0("MULTISCHOLAR_GUI_ENTRY_ORACLE_DIR=", shQuote(file.path(out, "oracles")))))

# (c) Comparison. Every side in the oracle's JSON form (one round trip, as guiEntryComparable makes it).
asOracleForm <- function(steps) jsonlite::fromJSON(jsonlite::toJSON(steps, auto_unbox = TRUE, null = "null", na = "string",
    digits = NA), simplifyVector = FALSE)
readSteps <- function(path) if (file.exists(path)) asOracleForm(jsonlite::read_json(path, simplifyVector = FALSE)$steps) else NULL
truncate <- function(x, n = 160L) {
    x <- gsub("[\t\r\n]+", " ", paste(as.character(unlist(x)), collapse = " "))
    if (nchar(x) > n) paste0(substr(x, 1L, n), "...") else x
}
isDigest <- function(value) is.character(value) && length(value) == 1L &&
    grepl("^((<table>|<widget>|<text>|<ui>) )?[0-9a-f]{64}$", value)
# One row per differing key (compareLane, packaging/verify/run-lanes.R:104-135); `a` is the first side, `b` the second.
compareSteps <- function(lane, first, second) {
    rows <- list()
    for (step in union(names(second), names(first))) {
        ours <- first[[step]]
        theirs <- second[[step]]
        if (is.null(ours) || is.null(theirs)) {
            rows[[length(rows) + 1L]] <- data.frame(lane = lane, step = step, field = "<step>", key = "", kind = "discrete",
                a = if (is.null(ours)) "<absent>" else "<present>", b = if (is.null(theirs)) "<absent>" else "<present>",
                stringsAsFactors = FALSE)
            next
        }
        for (field in union(names(ours), names(theirs))) {
            x0 <- ours[[field]]
            y0 <- theirs[[field]]
            if (identical(x0, y0)) next
            keyed <- is.list(x0) && is.list(y0) && !is.null(names(x0)) && !is.null(names(y0))
            for (key in if (keyed) union(names(x0), names(y0)) else "") {
                x <- if (keyed) x0[[key]] else x0
                y <- if (keyed) y0[[key]] else y0
                if (identical(x, y)) next
                rows[[length(rows) + 1L]] <- data.frame(lane = lane, step = step, field = field, key = truncate(key, 200L),
                    kind = if (isDigest(x) && isDigest(y)) "digest" else "discrete",
                    a = if (is.null(x)) "<absent>" else truncate(x), b = if (is.null(y)) "<absent>" else truncate(y),
                    stringsAsFactors = FALSE)
            }
        }
    }
    if (!length(rows)) return(NULL)
    do.call(rbind, rows)
}
classifyFrozenBase <- function(rows) {
    integration <- rows$field == "artifacts" & ((startsWith(rows$key, "integration/proteomics/") &
        (rows$a == "<absent>" | rows$b == "<absent>")) | rows$key == "<store entries>")
    ifelse(integration, "integration_artifacts", "outside")
}
classifyBaseHead <- function(rows, head_steps) {
    routing <- rows$step %in% routing_steps & (
        (rows$field == "tab_status" & rows$key %in% c("quality_control", "normalization")) |
        (rows$field == "outputs" & rows$key == "proteomics_workflow-workflow_progress"))
    ifelse(routing, "qc_routing (KIT-180)", "outside")
}

all_rows <- list()
step_rows <- list()
verdict <- list()
check <- function(name, value) verdict[[length(verdict) + 1L]] <<- data.frame(check = name, value = as.character(value),
    stringsAsFactors = FALSE)
reference <- utils::read.delim(file.path(head_root, "docs/qa/workflow-execution/kit/kit-180/g1a-reference-protein-lanes.tsv"),
    stringsAsFactors = FALSE)
stepDigest <- function(step) if (is.null(step)) "<absent>" else digest::digest(step, algo = "sha256")
for (lane in lanes) {
    frozen <- readSteps(file.path(base_tree, "tests", "testdata", "gui-entry", paste0(lane, ".json")))
    base <- readSteps(file.path(out, "base", paste0(lane, ".json")))
    head <- readSteps(file.path(out, "head", paste0(lane, ".json")))
    oracle <- readSteps(file.path(out, "oracles", paste0(lane, ".json")))
    check(paste(lane, "files_present"), !is.null(frozen) && !is.null(base) && !is.null(head) && !is.null(oracle))
    if (is.null(frozen) || is.null(base) || is.null(head) || is.null(oracle)) next
    # Step 1: base against the frozen file, G1(a) chain link (i).
    step1 <- compareSteps(lane, frozen, base)
    if (!is.null(step1)) {
        step1$class <- classifyFrozenBase(step1)
        step1$comparison <- "frozen_vs_base"
        all_rows[[length(all_rows) + 1L]] <- step1
    }
    classed <- if (is.null(step1)) character() else unique(step1$step[step1$class == "integration_artifacts"])
    outside1 <- if (is.null(step1)) character() else unique(step1$step[step1$class == "outside"])
    expected <- reference$step[reference$lane == lane & reference$class == "integration_artifacts"]
    check(paste(lane, "step1_outside_steps"), paste(outside1, collapse = ","))
    check(paste(lane, "step1_integration_artifacts_steps"), paste(sort(setdiff(classed, outside1)), collapse = ","))
    check(paste(lane, "step1_matches_reference"), length(outside1) == 0L && setequal(classed, expected))
    # Step 2: head against base, live; the changed step is `design`.
    step2 <- compareSteps(lane, base, head)
    if (!is.null(step2)) {
        step2$class <- classifyBaseHead(step2)
        step2$comparison <- "base_vs_head"
        all_rows[[length(all_rows) + 1L]] <- step2
    }
    routed <- vapply(routing_steps, function(step) identical(head[[step]]$tab_status$quality_control, "pending") &&
        identical(head[[step]]$tab_status$normalization, "disabled") &&
        identical(base[[step]]$tab_status$quality_control, "complete") &&
        identical(base[[step]]$tab_status$normalization, "pending"), logical(1))
    after <- head$qc_protein_replicate_filter$tab_status
    check(paste(lane, "step2_outside_rows"), if (is.null(step2)) 0L else sum(step2$class == "outside"))
    check(paste(lane, "step2_routing_rows"), if (is.null(step2)) 0L else sum(step2$class != "outside"))
    check(paste(lane, "step2_five_steps_routed"), all(routed))
    check(paste(lane, "step2_steps_differing"), if (is.null(step2)) "" else paste(unique(step2$step), collapse = ","))
    check(paste(lane, "t4_replicate_filter_complete_pending"), identical(after$quality_control, "complete") &&
        identical(after$normalization, "pending"))
    # The head record against the file frozen from it.
    step3 <- compareSteps(lane, head, oracle)
    check(paste(lane, "head_record_vs_oracle_rows"), if (is.null(step3)) 0L else nrow(step3))
    if (!is.null(step3)) {
        step3$class <- "outside"
        step3$comparison <- "head_vs_oracle"
        all_rows[[length(all_rows) + 1L]] <- step3
    }
    for (step in union(names(frozen), names(head))) {
        cause <- if (!is.null(step1) && step %in% step1$step) paste(unique(step1$class[step1$step == step]), collapse = ",") else ""
        step_rows[[length(step_rows) + 1L]] <- data.frame(lane = lane, step = step, frozen = stepDigest(frozen[[step]]),
            base = stepDigest(base[[step]]), head = stepDigest(head[[step]]), frozen_vs_base_cause = cause,
            routing_step = step %in% routing_steps, stringsAsFactors = FALSE)
    }
}
check("base_side_exit_status", base_status)
check("head_side_exit_status", head_status)
check("head_commit", head_commit)
check("base_commit", base_commit)
if (length(all_rows)) {
    table <- do.call(rbind, all_rows)
    table$changed_step <- "design"
    table$routing_step <- table$step %in% routing_steps
    table <- table[, c("comparison", "lane", "step", "field", "key", "kind", "a", "b", "class", "changed_step", "routing_step")]
    utils::write.table(table, file.path(out, "gui-entry-control.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
if (length(step_rows)) utils::write.table(do.call(rbind, step_rows), file.path(out, "gui-entry-steps.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE)
utils::write.table(do.call(rbind, verdict), file.path(out, "verdict.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(do.call(rbind, verdict), row.names = FALSE)
invisible(system2("git", c("worktree", "remove", "--force", shQuote(base_tree)), stdout = TRUE, stderr = TRUE))
quit(status = if (base_status == 0L && head_status == 0L) 0L else 1L)
