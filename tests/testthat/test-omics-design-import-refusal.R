# GUI audit R4, P-11 / A103-05 and A103-04: the design import's preflight, live for lipid and metab. An imported design
# the DA preflight finds invalid (no `group`, no `replicates`, a run the assays lack) is refused with
# "Invalid <omic> design import: ..." before the bus, the owner or the source directory is written. Before the repair
# the metab import ran no preflight: it wrote the configuration, mapping, design and contrasts, then failed.

designImportRefusalFixture <- function(omic) {
    if (identical(omic, "lipidomics")) {
        payload <- module_ci_lipid_design_payload(layout = "combined")
        list(payload = payload, workflow = module_ci_lipid_design_workflow(payload), paths = module_ci_lipid_design_paths(),
            write = module_ci_lipid_design_write_import_pack)
    } else {
        payload <- module_ci_metab_design_payload(layout = "combined")
        list(payload = payload, workflow = module_ci_metab_design_workflow(payload), paths = module_ci_metab_design_paths(),
            write = module_ci_metab_design_write_import_pack)
    }
}

# The bus fields the import writes, the owner's history and the files under the source directory.
designImportRefusalWritten <- function(fixture) {
    workflow <- fixture$workflow
    list(design_matrix = workflow$design_matrix, contrasts_tbl = workflow$contrasts_tbl, config_list = workflow$config_list,
        column_mapping = workflow$column_mapping, data_cln = getWorkflowOmicsBusTables(workflow, "data_cln"),
        data_tbl = getWorkflowOmicsBusTables(workflow, "data_tbl"), tab_status = workflow$tab_status,
        history = workflow$state_manager$getHistory(), pending = isDIASessionOwnerPending(workflow$state_manager),
        files = sort(list.files(fixture$paths$source_dir, recursive = TRUE, all.files = TRUE)))
}

# The import as the shared frame runs it (the reader, then the commit), through the real vroom reader.
designImportRefusalRun <- function(omic, fixture, edit = identity) {
    root <- withr::local_tempdir(pattern = "r4-design-import-", .local_envir = parent.frame())
    fixture$write(root, fixture$payload)
    edit(root)
    config <- validateOmicGuiConfig(omicGuiConfig(omic))
    config$design$design_import$qc_baseline <- function(...) invisible(NULL)
    result <- tryCatch(suppressWarnings(suppressMessages(shiny::withReactiveDomain(omicsOwnerSession(), {
        imported <- readOmicsDesignImport(root, config, fixture$paths$source_dir)
        commitOmicsDesignImport(fixture$workflow, imported, fixture$paths$source_dir, config)
    }))), error = function(error) error)
    list(result = result, root = root)
}

designImportDropColumn <- function(column) {
    function(root) {
        path <- file.path(root, "design_matrix.tab")
        design <- utils::read.delim(path, check.names = FALSE)
        utils::write.table(design[setdiff(names(design), column)], path, sep = "\t", quote = FALSE, row.names = FALSE)
    }
}

designImportRefusalEdits <- list(
    no_group = list(edit = designImportDropColumn("group"), error = "Design matrix is missing group column 'group'."),
    no_replicates = list(edit = designImportDropColumn("replicates"),
        error = "Design matrix is missing required 'replicates' column."),
    drifted_run = list(edit = function(root) {
        path <- file.path(root, "design_matrix.tab")
        design <- utils::read.delim(path, check.names = FALSE)
        design$Run[[1L]] <- "ZZ_9"
        utils::write.table(design, path, sep = "\t", quote = FALSE, row.names = FALSE)
    }, error = "ZZ_9"))

for (omic in c("lipidomics", "metabolomics")) local({
    omic <- omic
    test_that(sprintf("P-11: an invalid imported design is refused before anything is written (%s)", omic), {
        for (case in names(designImportRefusalEdits)) {
            fixture <- designImportRefusalFixture(omic)
            before <- designImportRefusalWritten(fixture)
            run <- designImportRefusalRun(omic, fixture, designImportRefusalEdits[[case]]$edit)
            expect_s3_class(run$result, "omics_design_refused")
            expect_match(conditionMessage(run$result), sprintf("Invalid %s design import: ", omic), fixed = TRUE, info = case)
            expect_match(conditionMessage(run$result), designImportRefusalEdits[[case]]$error, fixed = TRUE, info = case)
            expect_identical(designImportRefusalWritten(fixture), before, info = case)
            expect_identical(before$history, "initial", info = case)
            expect_true(before$pending, info = case)
        }
    })
})

test_that("P-11: a valid metab import is as before: the contrasts are the file's, not the preflight's normalised table", {
    fixture <- designImportRefusalFixture("metabolomics")
    run <- designImportRefusalRun("metabolomics", fixture)
    expect_false(inherits(run$result, "error"))
    expect_identical(fixture$workflow$contrasts_tbl, readOmicsDesignImportContrasts(file.path(run$root, "contrast_strings.tab")))
    expect_identical(fixture$workflow$state_manager$getHistory(), c("initial", "metab_raw_data_s4"))
})
