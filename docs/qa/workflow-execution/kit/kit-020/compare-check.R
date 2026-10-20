# KIT-020 audit F-01: join owner-check.json (kit_020, owner path) and control-check.json (kit_036, in memory).
`%||%` <- function(a, b) if (is.null(a)) b else a
dir <- "docs/qa/workflow-execution/kit/kit-020"
owner <- jsonlite::read_json(file.path(dir, "owner-check.json"))
control <- jsonlite::read_json(file.path(dir, "control-check.json"))
rows <- list()
for (name in c("maxquant", "fragpipe", "pd_tmt")) {
    o <- owner[[name]]; c <- control[[name]]
    for (state in names(c$state_hashes)) rows[[length(rows) + 1L]] <- data.frame(format = name, item = state,
        owner = o$state_hashes[[state]] %||% NA_character_, control = c$state_hashes[[state]], stringsAsFactors = FALSE)
    for (contrast in names(c$da$table_hashes)) for (table in names(c$da$table_hashes[[contrast]]))
        rows[[length(rows) + 1L]] <- data.frame(format = name, item = paste("DA", contrast, table),
            owner = o$da$table_hashes[[contrast]][[table]] %||% NA_character_, control = c$da$table_hashes[[contrast]][[table]], stringsAsFactors = FALSE)
}
table <- do.call(rbind, rows)
table$equal <- !is.na(table$owner) & table$owner == table$control
print(table[, c("format", "item", "equal")], row.names = FALSE)
cat(sprintf("%d of %d hashes equal; owner build %s (%s), control library %s (%s)\n", sum(table$equal), nrow(table),
    owner$head_commit, if (isTRUE(owner$tree_clean)) "clean tree" else "dirty tree", control$library$source_commit, control$library$path))
jsonlite::write_json(list(rows = table, equal = sum(table$equal), total = nrow(table), owner_commit = owner$head_commit,
    control_commit = control$library$source_commit), file.path(dir, "compare-check.json"), auto_unbox = TRUE, pretty = TRUE)
