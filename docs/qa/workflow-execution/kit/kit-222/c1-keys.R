# KIT-222 control C1 (ticket section 12): G1 (revised) and G6 on base and on the fix.
# Base is the base commit's technical-replicate edit (cf8af245, R/func_omics_design_shared.R:629-655) copied verbatim
# below (a pure function of the design; the metadata assignment is unchanged except for clearing `tech_reps`, which G1
# and G6 do not exercise); head is the package as installed or loaded from the working tree.
# Writes <out>/compare.tsv (sequence, mode, base keys, head keys, base merged, head merged, identical) for every G1
# sequence, and <out>/g6.tsv (side, states explored, depth, first violation).
# Usage: Rscript c1-keys.R <out directory>
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1L]] else "c1"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(".", export_all = TRUE, helpers = FALSE, quiet = TRUE)
} else {
    library(MultiScholaR)
}
ns <- asNamespace("MultiScholaR")
headEdit <- function(d, runs, mode, number = NULL) get("applyOmicsDesignTechRepAssignment", ns)(d, runs, mode, number)
assign_meta <- function(d, runs, group) get("applyOmicsDesignMetadataAssignment", ns)(d, runs, group, "", "", 1L)

baseEdit <- function(designMatrix, techRepSamples, assignmentMode, manualReplicateNumber = NULL) {
    refuse <- function(message) list(ok = FALSE, type = "warning", message = message, designMatrix = designMatrix)
    if (length(techRepSamples) < 2) return(refuse("Please select at least two samples."))
    selected <- which(designMatrix$Run %in% techRepSamples)
    selectedGroups <- unique(designMatrix$group[selected])
    if (length(selectedGroups) > 1 || any(is.na(selectedGroups))) {
        return(refuse("All selected samples must belong to the same group."))
    }
    baseReplicateNumber <- switch(assignmentMode,
        lowest = min(designMatrix$replicates[selected], na.rm = TRUE),
        first = designMatrix$replicates[designMatrix$Run == techRepSamples[1]],
        manual = manualReplicateNumber,
        stop("Unsupported technical replicate assignment mode."))
    originalReplicateNumbers <- unique(designMatrix$replicates[selected])
    consolidated <- length(originalReplicateNumbers) - 1
    highest <- max(originalReplicateNumbers, na.rm = TRUE)
    toAdjust <- which(designMatrix$group == selectedGroups[1] & !is.na(designMatrix$group) &
        designMatrix$replicates > highest & !designMatrix$Run %in% techRepSamples)
    if (length(toAdjust) > 0) designMatrix$replicates[toAdjust] <- designMatrix$replicates[toAdjust] - consolidated
    designMatrix$replicates[selected] <- baseReplicateNumber
    techRepMapping <- stats::setNames(seq_along(techRepSamples), techRepSamples)
    designMatrix$tech_reps[selected] <- techRepMapping[designMatrix$Run[selected]]
    list(ok = TRUE, type = "message", message = paste("Assigned", length(techRepSamples), "samples as technical replicates"),
        designMatrix = designMatrix)
}

fresh <- function(runs) {
    data.frame(Run = runs, group = NA_character_, factor1 = NA_character_, factor2 = NA_character_,
        factor3 = NA_character_, batch = NA_character_, replicates = NA_integer_, tech_reps = NA_integer_,
        stringsAsFactors = FALSE)
}
perms <- function(x) {
    if (length(x) <= 1L) return(list(x))
    out <- list()
    for (i in seq_along(x)) for (rest in perms(x[-i])) out[[length(out) + 1L]] <- c(x[i], rest)
    out
}
keysText <- function(d) paste(sprintf("%s:%s", paste(d$group, d$replicates, sep = "_"), d$Run), collapse = " ")
merged <- function(d) {
    sample <- sub("_[0-9]+$", "", d$Run)
    any(vapply(split(sample, paste(d$group, d$replicates)), function(s) length(unique(s)) > 1L, logical(1)))
}

rows <- list()
sequence <- 0L
record <- function(mode, base, head) {
    sequence <<- sequence + 1L
    rows[[length(rows) + 1L]] <<- data.frame(sequence = sequence, mode = mode, base_keys = keysText(base),
        head_keys = keysText(head), base_merged = merged(base), head_merged = merged(head), identical = identical(base, head))
}
run <- function(order_rows, samples, order, mode, typed = NULL) {
    base <- head <- assign_meta(fresh(order_rows), order_rows, "G")
    for (i in seq_along(order)) {
        selection <- order_rows[startsWith(order_rows, paste0(order[[i]], "_"))]
        if (identical(mode, "first")) selection <- selection[order(base$replicates[match(selection, base$Run)])]
        number <- if (identical(mode, "manual")) typed[[i]] else NULL
        base <- baseEdit(base, selection, mode, number)$designMatrix
        head <- headEdit(head, selection, mode, number)$designMatrix
    }
    record(mode, base, head)
}
for (k in 2:3) {
    samples <- LETTERS[seq_len(k)]
    runs <- as.vector(t(outer(samples, 1:2, paste, sep = "_")))
    for (r in perms(runs)) for (o in perms(samples)) {
        run(r, samples, o, "lowest")
        run(r, samples, o, "first")
        run(r, samples, o, "manual", seq_len(k))
    }
}
samples <- LETTERS[1:4]
runs <- as.vector(t(outer(samples, 1:2, paste, sep = "_")))
for (r in perms(runs)) run(r, samples, samples, "manual", 1:4)
samples <- sprintf("S%02d", 1:50)
runs <- as.vector(t(outer(samples, 1:2, paste, sep = "_")))
set.seed(222L)
for (i in seq_len(200L)) {
    r <- runs[sample.int(length(runs))]
    for (mode in c("lowest", "first", "manual")) run(r, samples, samples, mode, seq_along(samples))
}
compare <- do.call(rbind, rows)
utils::write.table(compare, file.path(out, "compare.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("G1: %d sequences, %d identical, base merged %d, head merged %d\n", nrow(compare), sum(compare$identical),
    sum(compare$base_merged), sum(compare$head_merged)))

# G6: one group of 5 runs numbered 1..5; every sequence of 3 edits (every selection of two or more runs, every mode,
# every typed number 1..5, every first-clicked run). Invariant after every successful edit: within a key, all marked
# runs carry the same marking action. States are deduplicated (design numbers, marks and the actions' partition).
g6 <- function(edit) {
    runs <- paste0("r", 1:5)
    start <- assign_meta(fresh(runs), runs, "G")
    selections <- unlist(lapply(2:5, function(k) utils::combn(runs, k, simplify = FALSE)), recursive = FALSE)
    edits <- list()
    for (s in selections) {
        edits[[length(edits) + 1L]] <- list(runs = s, mode = "lowest", number = NULL)
        for (f in s) edits[[length(edits) + 1L]] <- list(runs = c(f, setdiff(s, f)), mode = "first", number = NULL)
        for (n in 1:5) edits[[length(edits) + 1L]] <- list(runs = s, mode = "manual", number = n)
    }
    canon <- function(d, action) paste(d$replicates, d$tech_reps, match(action, unique(action)), collapse = "|")
    violates <- function(d, action) {
        marked <- !is.na(d$tech_reps)
        any(vapply(split(action[marked], d$replicates[marked]), function(a) length(unique(a)) > 1L, logical(1)))
    }
    frontier <- list(list(d = start, action = rep(NA_integer_, 5L), path = character()))
    seen <- new.env(hash = TRUE)
    explored <- 0L
    for (depth in 1:3) {
        nxt <- list()
        for (state in frontier) for (e in edits) {
            result <- edit(state$d, e$runs, e$mode, e$number)
            explored <- explored + 1L
            if (!isTRUE(result$ok)) next
            action <- state$action
            action[match(e$runs, state$d$Run)] <- depth * 1000L + length(nxt)
            path <- c(state$path, sprintf("%s(%s%s)", e$mode, paste(e$runs, collapse = ","),
                if (is.null(e$number)) "" else paste0(";", e$number)))
            if (violates(result$designMatrix, action)) {
                return(data.frame(explored = explored, depth = depth, violation = paste(path, collapse = " then ")))
            }
            key <- canon(result$designMatrix, action)
            if (!is.null(seen[[key]])) next
            seen[[key]] <- TRUE
            nxt[[length(nxt) + 1L]] <- list(d = result$designMatrix, action = action, path = path)
        }
        frontier <- nxt
    }
    data.frame(explored = explored, depth = 3L, violation = "none")
}
g6rows <- rbind(cbind(side = "base", g6(baseEdit)), cbind(side = "head", g6(headEdit)))
utils::write.table(g6rows, file.path(out, "g6.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(g6rows)
