# KIT-113 W1 probe P1 (R2), the tables (scratch, not frozen). Run on the start commit's R code:
#   Rscript docs/qa/workflow-execution/kit/kit-113/w1/probe-p1-tables.R <package source>
# The annotation tables and the organism mappings the request half would send, each through the two transfer codecs
# D3 names: a KIT-019 saved display table (`saveDIADisplayTable` -> `materializeWorkflowTable`) and an `rds` state value
# (`saveWorkflowStateValue(..., codec = "rds")` -> `resolveWorkflowStateValue`), on an on-disk DIA owner.

args <- commandArgs(TRUE)
tree <- normalizePath(args[[1L]])
setwd(tree)
suppressMessages(pkgload::load_all(tree, quiet = TRUE, export_all = TRUE, helpers = TRUE))
source(file.path(tree, "docs", "qa", "workflow-execution", "kit", "kit-113", "w1", "probe-lib.R"))

root <- normalizePath(tempfile("kit113-p1-"), mustWork = FALSE)
dir.create(root)
manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "kit113-p1", "DIA"), WorkflowState$new(audit_enabled = FALSE))
manager$.context$reference_states <- TRUE

design <- module_ci_prot_da_design()
rownames(design) <- NULL
quant <- module_ci_prot_da_table()
rownames(quant) <- NULL
object <- module_ci_prot_da_object(design = design, protein_quant_table = quant)
deterministic <- buildProtEnrichDeterministicUniprotAnnotations(object, "9606")
lane <- readRDS(file.path(tree, "tests", "testdata", "gui-entry", "uniprot", "uniprot_annotations.RDS"))
mapped <- resolveProtEnrichOrganismMapping(workflowData = list(mixed_species_analysis = list(enabled = TRUE,
    selected_organism = "Homo sapiens", organism_mapping = NULL)), uniprotDatCln = deterministic, targetTaxon = "9606",
    catFn = function(...) invisible(NULL))$organismMapping

tables <- list(
    `the lane's UniProt cache (tests/testdata/gui-entry/uniprot/uniprot_annotations.RDS)` = lane,
    `the deterministic table on W1 (b)'s object (buildProtEnrichDeterministicUniprotAnnotations)` = deterministic,
    `the deterministic table, default identifiers` = buildProtEnrichDeterministicUniprotAnnotations(NULL, "9606"),
    `W1 (b)'s organism_filter no_mapping fixture (zero rows)` = data.frame(uniprot_acc = character(), taxon_id = character()),
    `W1 (b)'s organism_filter mapped session: the mapping from the annotations` = mapped)

describe <- function(x) {
    if (is.null(x)) return("NULL")
    sprintf("class %s; %d x %d; types %s; rownames %s; attributes %s", paste(class(x), collapse = "/"), nrow(x), ncol(x),
        paste(vapply(x, function(column) paste(class(column), collapse = "/"), ""), collapse = ","),
        if (is.integer(attr(x, "row.names"))) "automatic" else "character",
        paste(setdiff(names(attributes(x)), c("names", "row.names", "class")), collapse = ",") %|% "none")
}
`%|%` <- function(a, b) if (!nzchar(a)) b else a

withDIAStateScope(manager, function() {
    for (name in names(tables)) {
        value <- tables[[name]]
        say("")
        say("P1 table: %s", name)
        say("  sent:          %s", describe(value))
        for (codec in c("table", "rds")) {
            reference <- tryCatch(kit113Save(manager, value, codec), error = function(e) e)
            if (inherits(reference, "error")) {
                say("  %-5s save refused: %s", codec, conditionMessage(reference))
                next
            }
            back <- if (is.null(reference)) NULL else kit113Load(manager, reference)
            say("  %-5s received: %s", codec, describe(back))
            say("  %-5s identical(): %s%s", codec, identical(value, back),
                if (!identical(value, back) && !is.null(back)) paste0("; all.equal: ",
                    paste(all.equal(value, back, check.attributes = TRUE), collapse = " / ")) else "")
            if (!is.null(reference)) say("  %-5s reference %s bytes (object.size), codec %s", codec,
                format(as.numeric(object.size(reference))), reference$payload$codec %||% "-")
        }
    }
})
closeDIAStateAdapter(manager)
unlink(root, recursive = TRUE)
