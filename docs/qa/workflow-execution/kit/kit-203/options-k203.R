# KIT-203 S13: every figure the option notes, the report sentences and NEWS show, regenerated from KIT-189's committed
# tables (results/*.tsv; score run 45311934798) and, for the figures computed after scoring, its committed rows
# (rows/*.tsv.gz), formatted as score-k189.R formats them (sprintf("%.1f", 100 * x)), and found in the shipped texts.
# Evidence code under docs/qa: no package code. Runs on a GitHub-hosted runner (run-tests `scripts`) with MultiScholaR
# installed from the bound commit, from the package root.
#
# Usage: Rscript docs/qa/workflow-execution/kit/kit-203/options-k203.R <out dir>
# Writes <out dir>/s13.tsv: per figure, what it is, the regenerated value, the text that must hold it, and whether it
# does; and s13-recorded.tsv: the figures recorded for the owner (not shown). Exit status 0 (a mismatch is a result).

args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
out <- args[[1L]]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
k189 <- file.path("docs", "qa", "workflow-execution", "kit", "kit-189")
read <- function(f) utils::read.delim(file.path(k189, "results", f), stringsAsFactors = FALSE)
pc <- function(x) sprintf("%.1f", 100 * x)
cl <- read("classifications.tsv")
ce <- read("cells.tsv")
pa <- read("paired.tsv")
deciding <- c(10L, 11L, 12L, 15L, 16L, 17L, 20L, 21L, 22L)
rules <- get("ENRICHMENT_SIGNIFICANCE_RULES", asNamespace("MultiScholaR"))
text <- function(column, rule) rules[[column]][rules$rule == rule]
news <- paste(readLines("NEWS.md", warn = FALSE, encoding = "UTF-8"), collapse = " ")
news <- gsub("\\s+", " ", news)

any_of <- function(arm) cl[cl$arm == arm & cl$estimand == "a_any_call", , drop = FALSE]
fdp_of <- function(arm) cl[cl$arm == arm & cl$estimand == "d_fdp", , drop = FALSE]
range_of <- function(x) paste0(pc(min(x)), " % to ", pc(max(x)), " %")
q_any <- any_of("Q")
q_fdp <- fdp_of("Q")
power <- function(arm) ce$power_planted[ce$arm == arm & ce$estimand == "d_fdp" & ce$cell %in% deciding]

# The rows: pooled shares over all significant terms of a setting, the paired planted calls
files <- list.files(file.path(k189, "rows"), pattern = "^s-.*\\.tsv\\.gz$", full.names = TRUE)
rows <- do.call(rbind, lapply(files, function(f) {
    x <- utils::read.delim(gzfile(f), colClasses = "character", na.strings = "NA", quote = "")
    x[as.integer(x$cell) >= 10L & x$status == "ok" & x$arm %in% c("Q", "BHf"),
        c("cell", "propagation", "replicate", "arm", "calls", "calls_null", "planted_tested", "planted_called")]
}))
rows$cell <- as.integer(rows$cell)
pooled <- stats::aggregate(cbind(null = as.numeric(calls_null), calls = as.numeric(calls)) ~ cell + propagation + arm,
    data = rows, FUN = sum)
pooled$share <- pooled$null / pooled$calls
pooled_max <- function(arm) max(pooled$share[pooled$arm == arm & pooled$cell %in% deciding & pooled$calls > 0])
tested <- rows[rows$planted_tested == "TRUE" & rows$cell %in% deciding, , drop = FALSE]
wide <- stats::reshape(tested[, c("cell", "propagation", "replicate", "arm", "planted_called")], direction = "wide",
    idvar = c("cell", "propagation", "replicate"), timevar = "arm")
wide$diff <- (wide$planted_called.Q == "TRUE") - (wide$planted_called.BHf == "TRUE")
paired_max <- function(prop) max(stats::aggregate(diff ~ cell, data = wide[wide$propagation == prop, ], FUN = mean)$diff)
bhf_only <- sum(wide$diff < 0)
cell13 <- ce$estimate[ce$cell == 13L & ce$propagation == "is_a_part_of_regulates" & ce$arm == "Q" &
    ce$estimand == "d_fdp"]
q_fails_prop <- q_any[q_any$classification == "fails" & q_any$propagation == "is_a_part_of_regulates", ]
q_fails_prop <- q_fails_prop[order(q_fails_prop$cell), ]

checks <- list(
    list("default note: random lists with any call", range_of(any_of("BHf")$estimate), "note", "bh_whole_family"),
    list("default note: planted mean FDP", range_of(fdp_of("BHf")$estimate), "note", "bh_whole_family"),
    list("default note: pooled share, worst setting (post hoc)", paste0("at most ", pc(pooled_max("BHf")), " %"), "note",
        "bh_whole_family"),
    list("default note: planted term called where tested", range_of(power("BHf")), "note", "bh_whole_family"),
    list("default note: paired, with parent terms (post hoc)",
        paste0("up to ", pc(paired_max("is_a_part_of_regulates")), " percentage points less"), "note", "bh_whole_family"),
    list("default note: paired, without parent terms (post hoc)", paste0("up to ", pc(paired_max("none")),
        " points less often without"), "note", "bh_whole_family"),
    list("default note: never more often (lists where only the default called the planted term)",
        if (bhf_only == 0L) "never more often" else sprintf("%d lists", bhf_only), "note", "bh_whole_family"),
    list("q-value note: random lists with any call", range_of(q_any$estimate), "note", "qvalue"),
    list("q-value note: settings failing 7.5 %", sprintf("in %d of 18 settings and above it without certainty in %d more",
        sum(q_any$classification == "fails"), sum(q_any$estimate > 0.075 & q_any$classification != "fails")), "note",
        "qvalue"),
    list("q-value note: failing with parent terms", sprintf("(E. coli full proteome %s %%, its 2,000-protein subset %s %%)",
        pc(q_fails_prop$estimate[[1L]]), pc(q_fails_prop$estimate[[2L]])), "note", "qvalue"),
    list("q-value note: 10 or more calls", paste0("Up to ", pc(max(ce$large_rate[ce$arm == "Q" &
        ce$estimand == "a_any_call"])), " %"), "note", "qvalue"),
    list("q-value note: planted mean FDP", sprintf("was up to %s %%, clearly above 7.5 %% in %d of 18 settings, all %s",
        pc(max(q_fdp$estimate)), sum(q_fdp$classification == "fails"),
        if (all(q_fdp$propagation[q_fdp$classification == "fails"] == "none")) "without parent terms" else "MIXED"),
        "note", "qvalue"),
    list("q-value note: pooled share, worst setting (post hoc)", paste0("it was up to ", pc(pooled_max("Q")), " %"),
        "note", "qvalue"),
    list("q-value note: lists of 25 with parent terms (post hoc)", paste0("the average was ", pc(cell13), " %"), "note",
        "qvalue"),
    list("q-value note: paired (post hoc)", sprintf("up to %s percentage points more often than the default with parent terms and up to %s points more often without",
        pc(paired_max("is_a_part_of_regulates")), pc(paired_max("none"))), "note", "qvalue"),
    list("default report: random lists with any call", paste0(range_of(any_of("BHf")$estimate), " of random lists"),
        "methods", "bh_whole_family"),
    list("default report: planted mean FDP", paste0("at most ", pc(max(fdp_of("BHf")$estimate)), " %"), "methods",
        "bh_whole_family"),
    list("q-value report: random lists with any call", range_of(q_any$estimate), "methods", "qvalue"),
    list("q-value report: settings failing", sprintf("in %d of 18 settings, %s of them with parent terms, both with lists of %s",
        sum(q_any$classification == "fails"), c("zero", "one", "two", "three")[nrow(q_fails_prop) + 1L],
        paste(unique(q_fails_prop$L), collapse = ",")), "methods", "qvalue"),
    list("q-value report: BH over the kept terms", paste0(range_of(any_of("BH")$estimate), " for Benjamini-Hochberg"),
        "methods", "qvalue"),
    list("NEWS: q-value without and with parent terms", sprintf("up to %s %% of random protein lists without parent terms and up to %s %% with them",
        pc(max(q_any$estimate[q_any$propagation == "none"])),
        pc(max(q_any$estimate[q_any$propagation == "is_a_part_of_regulates"]))), "news", ""))
table <- do.call(rbind, lapply(checks, function(x) {
    target <- if (identical(x[[3L]], "news")) news else text(x[[3L]], x[[4L]])
    data.frame(figure = x[[1L]], expected = x[[2L]], text = paste(x[[3L]], x[[4L]]),
        found = grepl(x[[2L]], target, fixed = TRUE), stringsAsFactors = FALSE)
}))
utils::write.table(table, file.path(out, "s13.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
recorded <- data.frame(figure = c("random-list settings with a q-value estimate above 0.075",
    "their cells and parent terms", "planted lists tested without parent terms (deciding cells)",
    "q-value pooled share, lists of 25 with parent terms (cell 13)", "q-value mean FDP among lists with a call, worst",
    "default mean FDP among lists with a call, worst"),
    value = c(sum(q_any$estimate > 0.075),
        paste(sprintf("%d %s", q_any$cell[q_any$estimate > 0.075], q_any$propagation[q_any$estimate > 0.075]),
            collapse = "; "),
        paste(ce$planted_tested[ce$arm == "BHf" & ce$estimand == "d_fdp" & ce$cell %in% deciding &
            ce$propagation == "none"], collapse = ","),
        pc(pooled$share[pooled$cell == 13L & pooled$propagation == "is_a_part_of_regulates" & pooled$arm == "Q"]),
        pc(max(stats::aggregate(share ~ cell + propagation, FUN = mean, data = transform(rows[rows$arm == "Q" &
            rows$cell %in% deciding & as.numeric(rows$calls) > 0, ], share = as.numeric(calls_null) / as.numeric(calls)))$share)),
        pc(max(stats::aggregate(share ~ cell + propagation, FUN = mean, data = transform(rows[rows$arm == "BHf" &
            rows$cell %in% deciding & as.numeric(rows$calls) > 0, ], share = as.numeric(calls_null) / as.numeric(calls)))$share))),
    stringsAsFactors = FALSE)
utils::write.table(recorded, file.path(out, "s13-recorded.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(table)
print(recorded)
writeLines(sprintf("S13: %d of %d figures found as regenerated", sum(table$found), nrow(table)))
