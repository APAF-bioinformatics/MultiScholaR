# KIT-193 section 8.2: the digests of target 3's bound objects and of target 1's fifth (the helper), on the installed
# package, for bindings.tsv (second proofs audit Q03 item 7, Q10 (a)); and the shipped P at (320, 137) for ticket
# 6.1 (iii) (Q10 (b)). Functions: `digest::digest(deparse(body(fn)), algo = "sha256")`; constants:
# `digest::digest(deparse(x), algo = "sha256")`. Run by run-tests after the package is installed.
# Usage: Rscript docs/qa/workflow-execution/kit/kit-193/scripts/bindings-t3.R <out dir>
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1L]] else "kit193-out"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
ns <- asNamespace("MultiScholaR")
fns <- c("goAncestorPairs", "propagateGoAnnotations", "buildUniprotGoGeneSets",
    "removeRowsWithMissingValuesPercentHelper")
consts <- c("goPropagationRelations", "goAspectRoots")
digests <- c(
    vapply(fns, function(name) digest::digest(deparse(body(get(name, envir = ns))), algo = "sha256"), character(1L)),
    vapply(consts, function(name) digest::digest(deparse(get(name, envir = ns)), algo = "sha256"), character(1L)))
rows <- data.frame(name = c(fns, consts), digest = unname(digests), stringsAsFactors = FALSE)
utils::write.table(rows, file.path(out, "bindings-t3.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# The shipped cutoff of a one-group design of 320 runs at min_reps_per_group 137, as test-kit193-pairs-filter.R
# computes it (helper-kit193-filter.R: `kit193Cutoffs`), with the replay marker's option unset.
methods::setClass("Kit193Design", slots = c(design_matrix = "data.frame", args = "list"))
option <- get("workflowReplayRuleMarkers", envir = ns)$missing_value_comparison$option
options(stats::setNames(list(NULL), option))
design <- data.frame(Run = paste0("G1_", seq_len(320L)), group = "G1", stringsAsFactors = FALSE)
object <- methods::new("Kit193Design", design_matrix = design, args = list())
env <- new.env()
assign("config_list", list(), envir = env)
object <- suppressMessages(get("updateMissingValueParameters", envir = ns)(object, min_reps_per_group = 137L,
    min_groups = 1L, grouping_variable = "group", env = env))
cutoff <- object@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff
p320 <- data.frame(nstar = 320L, m = 137L, cutoff = format(cutoff, digits = 17L),
    P = format(round(1000 * cutoff), scientific = FALSE))
utils::write.table(p320, file.path(out, "p-320-137.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(R.version.string, paste("MultiScholaR", format(utils::packageVersion("MultiScholaR"))),
    paste("digest", format(utils::packageVersion("digest")))), file.path(out, "versions.txt"))
print(rows)
print(p320)
