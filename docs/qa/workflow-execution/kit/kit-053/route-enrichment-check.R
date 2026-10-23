# KIT-053: the enrichment declaration on the DIA_LIMPA route's own DA manifest (a copy of the accepted screen
# lifecycle member's project): the declaration runs on the manifest's source state, saves its product and records the
# reported service in the head (KIT-026), and the analysis record carries it. The kernel is the package's deterministic
# double (the synthetic corpus's protein ids map in no external service), as test-workflow-dia-products.R uses it.
# Usage: KIT053_INSTALLED=1 R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-053/route-enrichment-check.R <project copy> <label>
suppressMessages(library(MultiScholaR)); ns <- asNamespace("MultiScholaR")
for (f in ls(ns)) if (!exists(f, envir = globalenv(), inherits = TRUE) && is.function(get(f, ns))) assign(f, get(f, ns), envir = globalenv())
a <- commandArgs(TRUE); project <- a[[1]]; label <- a[[2]]
harness <- new.env(parent = globalenv()); sys.source("tools/profiling/workflow_kit.R", harness)
store <- harness$kitFindStore(project); root <- dirname(store)
id <- readSavedStoreIdentity(file.path(root, ".workflow-entities", "store.json"))
manager <- newDIAStateAdapter(newWorkflowEntityStore(root, id$project_id, id$workflow_id), reopen = TRUE)
out <- withDIAStateScope(manager, function() {
  manifest <- resolveDIAStateValue(manager$.context, manager$.context$snapshot$da_results)
  long <- materializeDIATable(manager$.context, manifest$results[[1]]$tables$da_proteins_long, NULL)
  long$uniprot_acc <- long[[manifest$protein_id_column]]
  contrast <- manifest$contrasts$contrasts[[1]]
  # The module's argument, built from the route's own DA table and design (the class the enrichment tab passes).
  design <- readDIAExecutionComponents(manager, manifest$source, fields = "design_matrix")@design_matrix
  da_results <- methods::new("da_results_for_enrichment", contrasts = tibble::tibble(contrasts = contrast),
    da_data = stats::setNames(list(long[long$comparison == unique(long$comparison)[[1]], , drop = FALSE]), contrast),
    design_matrix = as.data.frame(design))
  arguments <- list(da_results = da_results, taxon_id = 9606, up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05,
    pathway_dir = tempfile("pathway-"), go_annotations = NULL, exclude_iea = FALSE, protein_id_column = "uniprot_acc",
    contrast_names = manifest$contrasts$contrasts, correction_method = "gSCS")
  dir.create(arguments$pathway_dir)
  before <- length(manager$.context$snapshot$external_resources)
  declared <- runDIAReferenceEnrichment(list2env(list(state_manager = manager)), manifest, arguments,
    process_fn = runProtEnrichDeterministicProcessEnrichments)
  resources <- manager$.context$snapshot$external_resources
  record <- buildWorkflowAnalysisRecord(manager)
  list(workflow_type = manager$workflow_type, route = diaActiveRoute(manager), da_model = manifest$model$model,
    service = declared$service, product_saved = isSavedWorkflowReference(declared$reference),
    head_resources_added = length(resources) - before,
    head_resource = resources[[length(resources)]][c("service", "release", "release_kind")],
    record_resources = length(record$inputs$external_resources %||% record$external_resources %||% list()),
    record_valid = isTRUE(tryCatch({ validateWorkflowAnalysisRecord(record); TRUE }, error = function(e) FALSE)))
})
out <- c(list(schema = "kit.route_enrichment_check.v1", ticket = "KIT-053", label = label, date = format(Sys.Date()),
  head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE)), out)
jsonlite::write_json(out, file.path("docs/qa/workflow-execution/kit/kit-053", paste0("route-enrichment-check-", label, ".json")), auto_unbox = TRUE, pretty = TRUE, null = "null")
cat("route", out$route, "| model", out$da_model, "| service", out$service, "| product", out$product_saved, "| head resources +", out$head_resources_added,
  "| record resources", out$record_resources, "| record valid", out$record_valid, "\n")
