# KIT-092: the KIT-005 journey oracles are recomputed from control-1's store by the bindings generator at every bind
# (KIT-012), except the products the driver attached to each step, which stay as recorded in-run. This checks those
# attached products against the store under the v1.6 and the v1.7 canonical form. Read-only.
# Usage, from the repository root: R_LIBS=<kit_065 library> Rscript kit-092/kit005-check.R <v1.6 workflow_kit.R>
P <- "/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/kit/runtime/kit-journey"
invisible(loadNamespace("MultiScholaR"))
api <- function(name) get(name, asNamespace("MultiScholaR"))
cur <- new.env(); sys.source("tools/profiling/workflow_kit.R", cur)
old <- new.env(); sys.source(commandArgs(TRUE)[[1L]], old)
result <- readRDS(file.path(P, "control-1", "journey-result.rds"))
cat("product_ids present:", sum(vapply(result$populations, function(p) !is.null(p$product_ids), TRUE)), "of", length(result$populations), "\n")
store <- unique(vapply(Filter(function(p) !is.null(p$state), result$populations), `[[`, "", "store"))
for (label in c("v16", "v17")) {
  h <- if (label == "v16") old else cur
  t0 <- Sys.time(); fp <- h$kitFingerprintStore(store, api)
  attached <- cur$kitStepProducts(result$populations, fp$products)
  diff <- 0L; total <- 0L
  for (step in names(result$populations)) {
    p <- result$populations[[step]]; if (is.null(p$fingerprint)) next
    inrun <- setdiff(names(p$fingerprint$products), names(fp$states[[p$state]]$products))
    total <- total + length(inrun)
    if (!setequal(inrun, names(attached[[step]]))) { diff <- diff + 1L; cat(label, step, "in-run", length(inrun), "recomputed", length(attached[[step]]), "only in-run:", head(setdiff(inrun, names(attached[[step]])), 3), "\n") }
  }
  cat(label, "attached products in-run", total, "steps differing", diff, round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
}
