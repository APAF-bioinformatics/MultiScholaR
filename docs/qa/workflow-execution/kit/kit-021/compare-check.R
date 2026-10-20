# KIT-021 acceptance criterion 2: join da-tables-<corpus>.json (kit_021, the oracle run's saved DA tables on the owner
# path, which the counted member reproduces hash-identically) with control-check-<corpus>.json (kit_036, in memory).
#   Rscript docs/qa/workflow-execution/kit/kit-021/compare-check.R lfq_300
`%||%` <- function(a, b) if (is.null(a)) b else a
args <- commandArgs(TRUE)
corpus_key <- if (length(args)) args[[1L]] else "lfq_300"
dir <- "docs/qa/workflow-execution/kit/kit-021"
owner <- jsonlite::read_json(file.path(dir, paste0("da-tables-", corpus_key, ".json")))
control <- jsonlite::read_json(file.path(dir, paste0("control-check-", corpus_key, ".json")))
stopifnot(identical(owner$corpus_sha256, control$corpus_sha256))
# The owner names a contrast by its DA-tab string (`groupA-groupB`), the control by the same string.
rows <- list()
for (contrast in names(control$da$table_hashes)) {
    o <- owner$table_hashes[[contrast]]
    for (table in names(control$da$table_hashes[[contrast]]))
        rows[[length(rows) + 1L]] <- data.frame(contrast = contrast, table = table,
            owner = o[[table]] %||% NA_character_, control = control$da$table_hashes[[contrast]][[table]], stringsAsFactors = FALSE)
}
table <- do.call(rbind, rows)
table$equal <- !is.na(table$owner) & table$owner == table$control
print(table[, c("contrast", "table", "equal")], row.names = FALSE)
missing <- setdiff(names(owner$table_hashes), names(control$da$table_hashes))
cat(sprintf("%d of %d DA table hashes equal across %d contrasts (%s owner contrasts without a control counterpart); owner build %s on %s, control %s on %s\n",
    sum(table$equal), nrow(table), length(unique(table$contrast)), length(missing), owner$build, owner$library$path, control$library$source_commit, control$library$path))
jsonlite::write_json(list(corpus = corpus_key, corpus_sha256 = owner$corpus_sha256, rows = table, equal = sum(table$equal), total = nrow(table),
    owner_rows = owner$rows, control_rows = control$da$rows, owner_build = owner$build, control_commit = control$library$source_commit,
    control_states = control$state_hashes, control_ruv = control$ruv, control_seconds = control$seconds),
    file.path(dir, paste0("compare-check-", corpus_key, ".json")), auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
