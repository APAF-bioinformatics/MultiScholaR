#!/usr/bin/env Rscript
# The KIT-080 rules on KIT-084's journey pairs, reusing the single-member receipts as control and candidate.
# Usage: Rscript pair-decisions.R <pair ticket label> <control stage directory> <candidate stage directory>
args <- commandArgs(TRUE)
repo <- "/home/doktersmol/Projects/APAF/MultiScholaR"
helpers <- new.env(parent = globalenv())
`%||%` <- function(a, b) if (is.null(a)) b else a
assign("%||%", `%||%`, envir = helpers)
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), envir = helpers)
bindings <- jsonlite::fromJSON(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"), simplifyVector = FALSE)
label <- args[[1L]]
control <- readRDS(file.path(args[[2L]], "stage.rds"))
candidate <- readRDS(file.path(args[[3L]], "stage.rds"))
bindings$tickets[[label]] <- list(rule = "whole_slice_anon", worker_steps = list("normalisation", "correlation", "da"), ceiling_seconds = 14400L)
receipt <- list(cases = list(`control-1` = control$cases[["candidate-1"]], `candidate-1` = candidate$cases[["candidate-1"]]))
selector <- paste0("kit:journey-", label)
cases <- list(list(id = "control-1", role = "control", selector = paste0(selector, ":control:1")),
    list(id = "candidate-1", role = "candidate", selector = paste0(selector, ":candidate:1")))
decision <- helpers$kitJourneyTicketDecision(receipt, cases, bindings)
decision$members <- NULL
output <- file.path(repo, "docs/qa/workflow-execution/kit/kit-084", paste0("pair-decision-", label, ".json"))
writeLines(jsonlite::toJSON(c(list(ticket = "KIT-084", pair = label, control = control$stage, candidate = candidate$stage), decision),
    auto_unbox = TRUE, pretty = TRUE, null = "null", digits = 6), output)
cat(label, decision$status, "failing", paste(decision$failing_rule %||% "none", collapse = ","), "\n")
