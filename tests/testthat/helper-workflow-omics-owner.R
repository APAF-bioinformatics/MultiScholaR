# KIT-070/071: the owner fixtures shared by test-workflow-omics-owner.R and test-workflow-omics-sessions.R: a fake
# Shiny session, a four-sample design and one small object per omic class.

omicsOwnerSession <- function() {
  session <- new.env(parent = emptyenv())
  session$userData <- new.env(parent = emptyenv())
  session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
  session
}

omicsOwnerDesign <- function() {
  data.frame(Run = c("S1", "S2", "S3", "S4"), group = c("A", "A", "B", "B"), stringsAsFactors = FALSE)
}

omicsOwnerLipidObject <- function(n = 6L) {
  assay <- function(prefix) {
    data.frame(lipid_id = paste0(prefix, "_L", seq_len(n)), lipid_name = paste("Lipid", seq_len(n)),
      S1 = seq_len(n) + 0.5, S2 = seq_len(n) + 1.5, S3 = seq_len(n) + 2.5, S4 = seq_len(n) + 3.5,
      check.names = FALSE, stringsAsFactors = FALSE)
  }
  methods::new("LipidomicsAssayData", lipid_data = list(pos = assay("pos"), neg = assay("neg")),
    design_matrix = omicsOwnerDesign(), sample_id = "Run", group_id = "group",
    lipid_id_column = "lipid_id", annotation_id_column = "lipid_name")
}

omicsOwnerMetabObject <- function(n = 5L) {
  assay <- function(prefix) {
    data.frame(metabolite_id = paste0(prefix, "_M", seq_len(n)), metabolite_name = paste("Met", seq_len(n)),
      S1 = seq_len(n) * 2, S2 = seq_len(n) * 3, S3 = seq_len(n) * 4, S4 = seq_len(n) * 5,
      check.names = FALSE, stringsAsFactors = FALSE)
  }
  methods::new("MetaboliteAssayData", metabolite_data = list(LCMS_Pos = assay("pos"), GCMS = assay("gc")),
    design_matrix = omicsOwnerDesign(), sample_id = "Run", group_id = "group",
    metabolite_id_column = "metabolite_id", annotation_id_column = "metabolite_name")
}


#' KIT-072: a bound lipidomics or metabolomics session owner in a temporary project for the module tests (the modules'
#' only state manager), with `object` committed under `state_name` when given. The owner is closed when the calling
#' test ends.
#' @param object A `LipidomicsAssayData` or `MetaboliteAssayData` to commit, or NULL for an empty bound owner.
#' @param state_name The state the object is committed under.
#' @param workflow The module bus the owner is attached to (an environment); a fresh one when NULL.
#' @param experiment_paths Paths for the owner; a temporary `source_dir` is added when missing.
#' @param bind FALSE leaves the owner pending, as at session start before the design step binds it.
#' @param env The frame whose exit closes the owner.
#' @return The bus, with `$state_manager` the owner.
omicsTestOwner <- function(object = NULL, state_name = NULL, workflow = NULL, experiment_paths = list(), config = list(),
    omic_type = if (methods::is(object, "MetaboliteAssayData")) "metabolomics" else "lipidomics", bind = TRUE,
    env = parent.frame()) {
  testthat::skip_if_not(workflowOwnerCapable())
  root <- withr::local_tempdir(.local_envir = env)
  if (is.null(experiment_paths$source_dir)) experiment_paths$source_dir <- file.path(root, "source")
  dir.create(experiment_paths$source_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(workflow)) workflow <- new.env(parent = emptyenv())
  if (is.null(workflow$config_list)) workflow$config_list <- config
  session <- omicsOwnerSession()
  workflow$state_manager <- createDIASessionOwner(workflow, experiment_paths, session, omic_type = omic_type)
  withr::defer(closeDIAStateAdapter(workflow$state_manager), envir = env)
  if (!isTRUE(bind)) return(workflow)
  shiny::withReactiveDomain(session, bindWorkflowOwnerTypeAtDesign(workflow,
    if (identical(omic_type, "metabolomics")) "METABOLOMICS" else "LIPIDOMICS"))
  if (!is.null(object)) workflow$state_manager$saveState(state_name %||% "test_state", object, workflow$config_list,
    "KIT-072 test fixture state")
  workflow
}

#' KIT-072: replace a declared wrapper (`runLipidReference*`, `runMetabReference*`, the store-step commits) for the
#' calling test, recording each call's arguments; module seam tests check the hand-off to the owner's one path while
#' the declarations' own tests check what the wrapper commits.
#' @param name The wrapper's name in the package namespace.
#' @param result The value the wrapper returns, or a function of the call's arguments.
#' @return An environment whose `calls` holds the recorded argument lists.
omicsMockWrapper <- function(name, result = NULL, env = parent.frame()) {
  capture <- new.env(parent = emptyenv())
  capture$calls <- list()
  mock <- function(...) {
    args <- list(...)
    capture$calls[[length(capture$calls) + 1L]] <- args
    if (is.function(result)) result(...) else result
  }
  do.call(testthat::local_mocked_bindings, c(stats::setNames(list(mock), name),
    list(.package = "MultiScholaR", .env = env)))
  capture
}
