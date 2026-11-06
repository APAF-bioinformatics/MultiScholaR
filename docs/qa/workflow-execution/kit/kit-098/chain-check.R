# KIT-098 README chain check, item 6 (the verifier's branch compares the member's state with its oracle), in process on the
# screen corpus: each recording's store is re-fingerprinted and compared with the bound oracles exactly as verifyKitMember's
# oracle branch does for a transcriptomics member (kitOracleJourneyPopulations, kitStepProducts, kitReproduceMismatches with
# the identity products, kitHistoryMismatches), after kitCases has accepted the case at binding. A positive control alters one
# oracle product and one history state and must be caught. Read-only (stores are opened, not written beyond their io log).
# Usage: R_LIBS=<kit_098 library> Rscript docs/qa/workflow-execution/kit/kit-098/chain-check.R
suppressMessages(library(MultiScholaR))
repo <- normalizePath(".")
ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
h <- new.env(parent = globalenv()); sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), h)
bindings <- h$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
for (stage in c("kit-journey-098a", "kit-journey-098b", "kit-lifecycle-trx", "kit-lifecycle-trx-fail")) {
    case <- h$kitCases(bindings, stage)[[1L]]
    cat(sprintf("kitCases %-24s accepted: route %s, oracles %s, in_memory %s, record %s\n", stage, case$route, case$oracles, case$in_memory, case$record))
}
for (ticket in c("098a", "098b")) {
    case <- h$kitCases(bindings, paste0("kit-journey-", ticket))[[1L]]
    recorded <- readRDS(file.path(repo, sprintf("docs/qa/workflow-execution/kit/kit-098/oracles-%s_trx_screen.rds", case$oracles)))
    fresh <- h$kitFingerprintStore(recorded$store, api)
    stopifnot(identical(fresh$head_sha256, recorded$head_sha256))
    attached <- h$kitStepProducts(recorded$populations, fresh$products, rendered = FALSE)
    populations <- recorded$populations
    for (step in names(populations)) {
        populations[[step]]$fingerprint <- fresh$states[[populations[[step]]$state]]
        populations[[step]]$fingerprint$products <- c(populations[[step]]$fingerprint$products, attached[[step]])
    }
    oracles <- h$kitOracleJourneyPopulations(bindings, case$route, case$oracles)
    mismatches <- c(h$kitReproduceMismatches(oracles, populations, ignore_products = h$kitOperationIdentityProducts),
        h$kitHistoryMismatches(bindings, case, h$kitReadHead(recorded$store, api)$history))
    judged <- names(Filter(function(o) !is.null(o$fingerprint), oracles))
    cat(sprintf("%s: steps judged %d (%s); mismatches %d%s\n", ticket, length(judged), paste(judged, collapse = ","), length(mismatches),
        if (length(mismatches)) paste0(" (", paste(mismatches, collapse = ", "), ")") else ""))
    # Positive control: one altered oracle product and one altered history state are caught.
    altered <- bindings
    key <- paste0(case$oracles, "_ruviii_export")
    name <- grep("^trx[.]ruviii:", names(altered$oracles[[key]]$products), value = TRUE)[[1L]]
    altered$oracles[[key]]$products[[name]]$sha256 <- "altered"
    altered$oracles[[paste0(case$oracles, "_history")]] <- rev(altered$oracles[[paste0(case$oracles, "_history")]])
    caught <- c(h$kitReproduceMismatches(h$kitOracleJourneyPopulations(altered, case$route, case$oracles), populations,
        ignore_products = h$kitOperationIdentityProducts), h$kitHistoryMismatches(altered, case, h$kitReadHead(recorded$store, api)$history))
    cat(sprintf("%s positive control: caught %s\n", ticket, paste(caught, collapse = ", ")))
    stopifnot(!length(mismatches), length(caught) == 2L)
}
cat("CHAIN ITEM 6 OK\n")
