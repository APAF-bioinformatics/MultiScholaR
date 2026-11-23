# KIT-177 script S2: the control and the freeze of the GUI-entry references prot_tmt, prot_lfq and prot_lfq_fragpipe
# (ticket section 6; KIT-159 G1(a) chain links (i) to (iii), ruling C75). In R0 the same script answers F4 on one lane
# (ticket section 8): there the head is the R0 commit, whose R/ and tests/ equal the base, so every comparison is
# expected to be empty.
#
# Usage, from the package root of the commit under test on a run-tests runner (full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-177/gui-entry-control.R <base commit> [<out dir>] [<lanes>]
#   <out dir>  default docs/qa/workflow-execution/kit/kit-177/control (made absolute)
#   <lanes>    comma-separated, default prot_tmt,prot_lfq,prot_lfq_fragpipe
#
# (a) Base side: a git worktree of <base commit>; a fresh Rscript process whose working directory is that tree's root,
#     under C collation from its start (as testthat sets it: withr::local_collate("C"); KIT-172 C5, KIT-180 R4b),
#     loads it (pkgload::load_all, helpers sourced into the global environment, no helper file edited) and runs
#     guiEntryRunLane() for each lane inside testthat::test_that() with MULTISCHOLAR_GUI_ENTRY_RECORD_DIR=<out>/base.
#     The helper guiEntryStateDigests is wrapped in the global environment: it returns the same digests and writes each
#     protein state's quantity table, once per state name and digest, to <out>/base-states/<lane>/ (all digits).
# (b) Head side: the same at this tree's root with guiEntryExpectLane() inside test_that(),
#     MULTISCHOLAR_GUI_ENTRY_RECORD_DIR=<out>/head, MULTISCHOLAR_GUI_ENTRY_FREEZE=1 and
#     MULTISCHOLAR_GUI_ENTRY_ORACLE_DIR=<out>/oracles; state tables to <out>/head-states/<lane>/.
# (c) Comparison, here: one row per differing key (packaging/verify/run-lanes.R compareLane's form):
#       frozen_vs_base  step 1 (link (i)): expected none, the files being KIT-180's re-freeze on run-tests; any row is
#                       `unexplained` (the ticket stops).
#       base_vs_head    steps 2 and 3 (links (ii), (iii)); the changed step is `design` on every lane. Before it every
#                       row is `unexplained`. From it on, a `states`, `artifacts` or `outputs` digest that differs on a key
#                       both sides hold is `plain_log2 (KIT-177)`; a key on one side only, or a difference in the head,
#                       history, tab status, notifications or refusal check (a step outcome, not a value) is
#                       `unexplained`.
#       head_vs_oracle  the head record against the file frozen from it (expected none).
#     State tables (step 3, cell by cell): every table both sides exported, by lane, state name and order of appearance.
#     At `protein_s4_initial` each cell is checked against the rule base = log2(v + 1), head = log2(v): the base minus the
#     head equals log2(1 + 2^-head) to 1e-12, or the base is 0 and the head missing (a zero turned missing); any other
#     cell is `unexplained`. For later tables the largest absolute difference and any change of shape or identifiers.
#     Written to <out>: gui-entry-control.tsv (comparison, lane, step, field, key, kind, a, b, class, changed_step,
#     before_changed_step), gui-entry-steps.tsv (lane, step, frozen/base/head step digests, before_changed_step),
#     gui-entry-tables.tsv (lane, state, occurrence, rows, cols, cells compared, largest absolute difference, design-rule
#     violations, zeros turned missing, other missing-value changes, class), verdict.tsv (check, value).
# Exit status non-zero when a side's process fails. Only <out>/oracles/<lane>.json for the three lanes are committed,
# and only when scripts.tsv shows this script `ok`, no `withheld.txt` exists and the verdict holds.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L, grepl("^[0-9a-f]{40}$", args[[1L]]))
base_commit <- args[[1L]]
out <- if (length(args) >= 2L) args[[2L]] else "docs/qa/workflow-execution/kit/kit-177/control"
lanes <- if (length(args) >= 3L) strsplit(args[[3L]], ",", fixed = TRUE)[[1L]] else c("prot_tmt", "prot_lfq", "prot_lfq_fragpipe")
changed_step <- "design"
for (sub in c("base", "head", "oracles", "base-states", "head-states")) dir.create(file.path(out, sub), recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
head_root <- normalizePath(".")
head_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
cat("KIT-177 S2: head", head_commit, "base", base_commit, "lanes", paste(lanes, collapse = ","), "out", out, "\n")
rscript <- file.path(R.home("bin"), "Rscript")

child <- file.path(tempdir(), "kit177-s2-child.R")
writeLines(c(
    "invisible(Sys.setlocale('LC_COLLATE', 'C'))",
    "args <- commandArgs(trailingOnly = TRUE)",
    "setwd(args[[1L]])",
    "side <- args[[2L]]",
    "lanes <- strsplit(args[[3L]], ',', fixed = TRUE)[[1L]]",
    "state_root <- args[[4L]]",
    "cat('KIT-177 S2 child', side, 'collate', Sys.getlocale('LC_COLLATE'), 'test_path',",
    "    normalizePath(file.path(getwd(), 'tests', 'testthat')), '\\n')",
    "suppressPackageStartupMessages(pkgload::load_all('.', quiet = TRUE, export_all = TRUE, helpers = FALSE))",
    "suppressMessages(testthat::source_test_helpers('tests/testthat', env = globalenv()))",
    "cat('KIT-177 S2 child', side, 'oracle dir', normalizePath(guiEntryOracleDir(), mustWork = FALSE), '\\n')",
    "exported <- new.env(parent = emptyenv())",
    "current_lane <- NA_character_",
    "original_digests <- guiEntryStateDigests",
    "guiEntryStateDigests <- function(manager, cache) {",
    "    digests <- original_digests(manager)",
    "    dir <- file.path(state_root, current_lane)",
    "    dir.create(dir, recursive = TRUE, showWarnings = FALSE)",
    "    for (name in names(digests)) {",
    "        key <- paste(current_lane, name, digests[[name]])",
    "        if (!is.null(exported[[key]])) next",
    "        exported[[key]] <- TRUE",
    "        object <- tryCatch(withWorkflowStateScope(manager, function() manager$getState(name)), error = function(e) NULL)",
    "        if (!isS4(object) || !'protein_quant_table' %in% methods::slotNames(object)) next",
    "        table <- as.data.frame(object@protein_quant_table)",
    "        numeric <- vapply(table, is.numeric, logical(1))",
    "        table[numeric] <- lapply(table[numeric], function(x) ifelse(is.na(x), 'NA', sprintf('%.17g', x)))",
    "        occurrence <- sum(startsWith(ls(exported), paste(current_lane, name, ''))) ",
    "        utils::write.table(table, file.path(dir, sprintf('%s--%d.tsv', name, occurrence)), sep = '\\t', quote = FALSE,",
    "            row.names = FALSE)",
    "    }",
    "    digests",
    "}",
    "failed <- FALSE",
    "for (lane in lanes) {",
    "    current_lane <- lane",
    "    ok <- tryCatch(if (identical(side, 'head')) {",
    "        testthat::test_that(paste('KIT-177 control head', lane), guiEntryExpectLane(lane))",
    "    } else {",
    "        testthat::test_that(paste('KIT-177 control base', lane), { guiEntryRunLane(lane); testthat::succeed() })",
    "    }, error = function(e) {",
    "        cat('lane', lane, 'error', class(e)[[1L]], conditionMessage(e), '\\n')",
    "        FALSE",
    "    })",
    "    cat('KIT-177 S2', side, lane, if (isTRUE(ok)) 'ok' else 'FAILED', '\\n')",
    "    if (!isTRUE(ok)) failed <- TRUE",
    "}",
    "quit(status = if (failed) 1L else 0L)"), child)
runSide <- function(tree, side, env) {
    console <- suppressWarnings(system2(rscript, c(shQuote(child), shQuote(tree), side, paste(lanes, collapse = ","),
        shQuote(file.path(out, paste0(side, "-states")))), stdout = TRUE, stderr = TRUE, env = env))
    cat(grep("^KIT-177 S2|^lane ", console, value = TRUE), sep = "\n")
    cat(utils::tail(console, 15L), sep = "\n")
    list(status = attr(console, "status") %||% 0L, console = console)
}

verdict <- list()
check <- function(name, value) verdict[[length(verdict) + 1L]] <<- data.frame(check = name,
    value = paste(as.character(value), collapse = ","), stringsAsFactors = FALSE)

# (a) Base side.
base_tree <- file.path(tempdir(), "kit177-s2-base")
invisible(system2("git", c("worktree", "add", "--detach", shQuote(base_tree), base_commit), stdout = TRUE, stderr = TRUE))
check("base_worktree_present", dir.exists(file.path(base_tree, "tests", "testthat")))
base_run <- runSide(base_tree, "base", paste0("MULTISCHOLAR_GUI_ENTRY_RECORD_DIR=", shQuote(file.path(out, "base"))))

# (b) Head side.
head_run <- runSide(head_root, "head", c(
    paste0("MULTISCHOLAR_GUI_ENTRY_RECORD_DIR=", shQuote(file.path(out, "head"))),
    "MULTISCHOLAR_GUI_ENTRY_FREEZE=1",
    paste0("MULTISCHOLAR_GUI_ENTRY_ORACLE_DIR=", shQuote(file.path(out, "oracles")))))
childLine <- function(run, pattern) sub(pattern, "", grep(paste0("^", pattern), run$console, value = TRUE))
check("base_child_collate", childLine(base_run, "KIT-177 S2 child base collate ")[1L])
check("head_child_collate", childLine(head_run, "KIT-177 S2 child head collate ")[1L])
check("base_test_path_inside_worktree", any(grepl(normalizePath(base_tree), childLine(base_run, "KIT-177 S2 child base collate "),
    fixed = TRUE)))
check("head_oracle_dir_redirected", any(grepl(file.path(out, "oracles"), childLine(head_run, "KIT-177 S2 child head oracle dir "),
    fixed = TRUE)))

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
compareSteps <- function(lane, first, second) {
    rows <- list()
    for (step in union(names(first), names(second))) {
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
stepIds <- function(steps) names(steps)
classifyBaseHead <- function(rows, order) {
    before <- match(rows$step, order) < match(changed_step, order)
    value_field <- rows$field %in% c("states", "artifacts", "outputs")
    both_sides <- rows$a != "<absent>" & rows$b != "<absent>"
    ifelse(!is.na(before) & !before & value_field & both_sides, "plain_log2 (KIT-177)", "unexplained")
}

# State tables, cell by cell.
readStateTable <- function(path) {
    table <- utils::read.delim(path, colClasses = "character", check.names = FALSE, na.strings = character())
    table
}
compareTables <- function(lane) {
    base_dir <- file.path(out, "base-states", lane)
    head_dir <- file.path(out, "head-states", lane)
    files <- intersect(list.files(base_dir, pattern = "\\.tsv$"), list.files(head_dir, pattern = "\\.tsv$"))
    only <- setdiff(union(list.files(base_dir, pattern = "\\.tsv$"), list.files(head_dir, pattern = "\\.tsv$")), files)
    rows <- list()
    for (file in only) rows[[length(rows) + 1L]] <- data.frame(lane = lane, state = sub("--[0-9]+\\.tsv$", "", file),
        occurrence = sub("^.*--([0-9]+)\\.tsv$", "\\1", file), rows = NA, cols = NA, cells = NA, max_abs_difference = NA,
        design_rule_violations = NA, zero_to_missing = NA, other_missing_changes = NA,
        class = "unexplained (exported on one side only)", stringsAsFactors = FALSE)
    for (file in files) {
        state <- sub("--[0-9]+\\.tsv$", "", file)
        base <- readStateTable(file.path(base_dir, file))
        head <- readStateTable(file.path(head_dir, file))
        same_shape <- identical(dim(base), dim(head)) && identical(names(base), names(head))
        id <- names(base)[[1L]]
        same_ids <- same_shape && identical(base[[id]], head[[id]])
        numeric <- names(base)[vapply(base, function(x) all(x == "NA" | !is.na(suppressWarnings(as.numeric(x)))), logical(1))]
        numeric <- setdiff(numeric, id)
        b <- if (same_ids) suppressWarnings(as.matrix(as.data.frame(lapply(base[numeric], as.numeric)))) else NULL
        h <- if (same_ids) suppressWarnings(as.matrix(as.data.frame(lapply(head[numeric], as.numeric)))) else NULL
        cells <- if (same_ids) length(b) else NA
        max_diff <- if (same_ids && any(!is.na(b) & !is.na(h))) max(abs(b - h), na.rm = TRUE) else if (same_ids) 0 else NA
        zero_to_missing <- if (same_ids) sum(!is.na(b) & b == 0 & is.na(h)) else NA
        other_missing <- if (same_ids) sum(is.na(b) != is.na(h)) - zero_to_missing else NA
        violations <- NA
        if (same_ids && identical(state, "protein_s4_initial")) {
            both <- !is.na(b) & !is.na(h)
            violations <- sum(abs((b[both] - h[both]) - log2(1 + 2^(-h[both]))) > 1e-12)
        }
        class <- if (!same_ids) "unexplained (shape or identifiers differ)" else
            if (identical(state, "protein_s4_initial") && (violations > 0L || other_missing != 0L)) "unexplained (design rule)" else
            if (max_diff == 0 && zero_to_missing == 0L && other_missing == 0L) "identical" else "plain_log2 (KIT-177)"
        rows[[length(rows) + 1L]] <- data.frame(lane = lane, state = state, occurrence = sub("^.*--([0-9]+)\\.tsv$", "\\1", file),
            rows = nrow(base), cols = ncol(base), cells = cells, max_abs_difference = signif(max_diff, 6L),
            design_rule_violations = violations, zero_to_missing = zero_to_missing, other_missing_changes = other_missing,
            class = class, stringsAsFactors = FALSE)
    }
    if (!length(rows)) return(NULL)
    do.call(rbind, rows)
}

all_rows <- list()
step_rows <- list()
table_rows <- list()
orders <- list()
stepDigest <- function(step) if (is.null(step)) "<absent>" else digest::digest(step, algo = "sha256")
for (lane in lanes) {
    frozen <- readSteps(file.path(base_tree, "tests", "testdata", "gui-entry", paste0(lane, ".json")))
    base <- readSteps(file.path(out, "base", paste0(lane, ".json")))
    head <- readSteps(file.path(out, "head", paste0(lane, ".json")))
    oracle <- readSteps(file.path(out, "oracles", paste0(lane, ".json")))
    check(paste(lane, "files_present"), !is.null(frozen) && !is.null(base) && !is.null(head) && !is.null(oracle))
    if (is.null(frozen) || is.null(base) || is.null(head) || is.null(oracle)) next
    order <- stepIds(frozen)
    orders[[lane]] <- order
    check(paste(lane, "step_lists_equal"), identical(stepIds(base), order) && identical(stepIds(head), order))
    step1 <- compareSteps(lane, frozen, base)
    if (!is.null(step1)) {
        step1$class <- "unexplained"
        step1$comparison <- "frozen_vs_base"
        all_rows[[length(all_rows) + 1L]] <- step1
    }
    check(paste(lane, "step1_rows"), if (is.null(step1)) 0L else nrow(step1))
    step23 <- compareSteps(lane, base, head)
    if (!is.null(step23)) {
        step23$class <- classifyBaseHead(step23, order)
        step23$comparison <- "base_vs_head"
        all_rows[[length(all_rows) + 1L]] <- step23
    }
    before <- if (is.null(step23)) 0L else sum(match(step23$step, order) < match(changed_step, order), na.rm = TRUE)
    check(paste(lane, "step2_rows_before_changed_step"), before)
    check(paste(lane, "step3_rows_plain_log2"), if (is.null(step23)) 0L else sum(step23$class == "plain_log2 (KIT-177)"))
    check(paste(lane, "step3_rows_unexplained"), if (is.null(step23)) 0L else sum(step23$class == "unexplained"))
    check(paste(lane, "base_vs_head_steps_differing"), if (is.null(step23)) "" else paste(unique(step23$step), collapse = ","))
    tables <- compareTables(lane)
    if (!is.null(tables)) table_rows[[length(table_rows) + 1L]] <- tables
    check(paste(lane, "state_tables_compared"), if (is.null(tables)) 0L else nrow(tables))
    check(paste(lane, "state_tables_unexplained"), if (is.null(tables)) 0L else sum(startsWith(tables$class, "unexplained")))
    check(paste(lane, "design_state_compared"), !is.null(tables) && "protein_s4_initial" %in% tables$state)
    step4 <- compareSteps(lane, head, oracle)
    check(paste(lane, "head_record_vs_oracle_rows"), if (is.null(step4)) 0L else nrow(step4))
    if (!is.null(step4)) {
        step4$class <- "unexplained"
        step4$comparison <- "head_vs_oracle"
        all_rows[[length(all_rows) + 1L]] <- step4
    }
    for (step in union(order, names(head))) step_rows[[length(step_rows) + 1L]] <- data.frame(lane = lane, step = step,
        frozen = stepDigest(frozen[[step]]), base = stepDigest(base[[step]]), head = stepDigest(head[[step]]),
        before_changed_step = isTRUE(match(step, order) < match(changed_step, order)), stringsAsFactors = FALSE)
}
check("checkout_oracles_unchanged", !length(system2("git", c("status", "--porcelain", "--", "tests/testdata/gui-entry"),
    stdout = TRUE)))
check("base_side_exit_status", base_run$status)
check("head_side_exit_status", head_run$status)
check("head_commit", head_commit)
check("base_commit", base_commit)
check("changed_step", changed_step)
if (length(all_rows)) {
    table <- do.call(rbind, all_rows)
    table$changed_step <- changed_step
    table$before_changed_step <- vapply(seq_len(nrow(table)), function(i) {
        order <- orders[[table$lane[[i]]]]
        isTRUE(match(table$step[[i]], order) < match(changed_step, order))
    }, logical(1))
    table <- table[, c("comparison", "lane", "step", "field", "key", "kind", "a", "b", "class", "changed_step", "before_changed_step")]
    utils::write.table(table, file.path(out, "gui-entry-control.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
if (length(step_rows)) utils::write.table(do.call(rbind, step_rows), file.path(out, "gui-entry-steps.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE)
if (length(table_rows)) utils::write.table(do.call(rbind, table_rows), file.path(out, "gui-entry-tables.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE)
utils::write.table(do.call(rbind, verdict), file.path(out, "verdict.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(do.call(rbind, verdict), row.names = FALSE)
invisible(system2("git", c("worktree", "remove", "--force", shQuote(base_tree)), stdout = TRUE, stderr = TRUE))
quit(status = if (base_run$status == 0L && head_run$status == 0L) 0L else 1L)
