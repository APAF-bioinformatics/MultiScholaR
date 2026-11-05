# KIT-073: the observable surface of the lipidomics and metabolomics declarations, dumped so that the shared core can
# be compared with the two files it replaced (tests/testthat/fixtures/kit073-omics-declarations-frozen.rds, taken on
# the base commit by docs/qa/workflow-execution/kit/kit-073/freeze-dump.R): every registry entry and declaration field
# (the kernel excepted), and both omics' wrapper journeys on the declaration fixtures, every committed state and
# product by content with times removed.

omicsDeclarationsScrub <- function(x) {
  if (inherits(x, c("POSIXt", "Date"))) return("<time>")
  if (is.function(x) || is.environment(x)) return("<fn>")
  if (isS4(x)) {
    for (slot in methods::slotNames(x)) methods::slot(x, slot, check = FALSE) <- omicsDeclarationsScrub(methods::slot(x, slot))
    return(x)
  }
  if (isSavedWorkflowReference(x)) return(list(reference = x$kind %||% "reference"))
  if (is.data.frame(x)) return(as.data.frame(lapply(x, omicsDeclarationsScrub), stringsAsFactors = FALSE, check.names = FALSE))
  if (is.list(x)) { x[] <- lapply(x, omicsDeclarationsScrub); return(x) }
  x
}

omicsDeclarationsRegistry <- function(omic_type) {
  names <- Filter(function(n) identical(workflowOperationEntry(n)$omic_type, omic_type), workflowOperationNames())
  stats::setNames(lapply(sort(names), function(n) {
    entry <- workflowOperationEntry(n)
    declaration <- entry$declaration
    entry$declaration <- NULL
    list(entry = entry, declaration = declaration[setdiff(names(declaration), "kernel")])
  }), sort(names))
}

omicsDeclarationsJourneyDump <- function(fixture, journey) {
  owner <- fixture$owner
  out <- journey(fixture)
  states <- owner$getHistory()
  list(
    history = states,
    states = stats::setNames(lapply(states, function(s) omicsDeclarationsScrub(owner$getState(s))), states),
    values = lapply(out[setdiff(names(out), "da")], function(r) omicsDeclarationsScrub(r$value)),
    correlations = omicsDeclarationsScrub(out$correlation$correlations),
    da = omicsDeclarationsScrub(list(manifest = out$da$manifest[setdiff(names(out$da$manifest), "source")],
      results = out$da$results[setdiff(names(out$da$results), "theObject")]))
  )
}

#' The dump of both omics: the registry and the journeys on fresh owners in `root`.
omicsDeclarationsDump <- function(root) {
  lipid <- lipidDeclarationsOwner(file.path(root, "lipid"))
  metab <- metabDeclarationsOwner(file.path(root, "metab"))
  on.exit({ closeDIAStateAdapter(lipid$owner); closeDIAStateAdapter(metab$owner) }, add = TRUE)
  list(
    registry = list(lipidomics = omicsDeclarationsRegistry("lipidomics"), metabolomics = omicsDeclarationsRegistry("metabolomics")),
    lipidomics = omicsDeclarationsJourneyDump(lipid, lipidDeclarationsJourney),
    metabolomics = omicsDeclarationsJourneyDump(metab, metabDeclarationsJourney)
  )
}

# KIT-089: the intensity filter defaults to `partitioned`, which hands the kernel only an empty template. A module test that
# fakes the kernel and checks what it received and returned characterises the `whole` path, so it pins that variant here.
kit089WholeIntensity <- function(env = parent.frame()) {
  original <- workflowOperationEntry
  testthat::local_mocked_bindings(workflowOperationEntry = function(operation) {
    entry <- original(operation)
    if (grepl("_intensity_filter$", operation)) entry$default_variant <- "whole"
    entry
  }, .package = "MultiScholaR", .env = env)
}
