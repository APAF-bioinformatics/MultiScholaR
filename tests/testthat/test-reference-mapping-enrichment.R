# KIT-167 section 10 and 8.3: proteomics enrichment with the reference identifiers.

k167Mapping <- data.frame(strain_accession = c("WP_1", "WP_2", "WP_3", "WP_4", "WP_5"),
    reference_accession_used = c("P1", "P1", "P3", "P4", NA),
    map_status = c("one_to_one", "one_to_one", "one_to_one", "one_to_one", "below_threshold"),
    stringsAsFactors = FALSE)

k167Da <- function(table) {
    methods::new("da_results_for_enrichment", contrasts = tibble::tibble(contrast = "B_vs_A"),
        da_data = list(B_vs_A = table), design_matrix = data.frame(sample = "S1", stringsAsFactors = FALSE))
}

test_that("the protein-group join: wrapper removed, human members ignored, ambiguous and unmapped groups", {
    groups <- c("ref|WP_1|;ref|WP_2|", "WP_3;sp|P69905|HBA_HUMAN", "P69905", "WP_1;WP_4", "WP_5", "ref|WP_3|")
    joined <- referenceAccessionForGroups(groups, k167Mapping)
    expect_identical(joined$reference_accession_used, c("P1", "P3", NA, NA, NA, "P3"))
    expect_identical(joined$group_class, c("mapped", "mapped", "unmapped", "ambiguous_group", "unmapped", "mapped"))
})

test_that("NA identifiers leave the query lists and the background; the 2 x 2 counts are per contrast", {
    table <- data.frame(Protein.IDs = c("WP_1", "WP_3", "WP_4;WP_1", "WP_5", "P69905"),
        fdr_qvalue = c(0.01, 0.5, 0.01, 0.01, 0.5), log2FC = c(2, 0, -2, 3, 0), stringsAsFactors = FALSE)
    attached <- attachReferenceMapping(k167Da(table), k167Mapping, "Protein.IDs", q_cutoff = 0.05)
    kept <- attached$da_results@da_data$B_vs_A
    expect_identical(kept$Protein.IDs, c("WP_1", "WP_3"))
    expect_false(anyNA(kept$reference_accession_used))
    counts <- attached$counts$B_vs_A
    expect_identical(c(counts$changed_mapped, counts$changed_unmapped, counts$unchanged_mapped,
        counts$unchanged_unmapped), c(1L, 2L, 1L, 1L))
    expect_identical(counts$ambiguous_group, 1L)
    expect_identical(counts$unmapped_changed_groups, c("WP_4;WP_1", "WP_5"))
    expect_identical(counts$mapped + counts$ambiguous_group + counts$unmapped, counts$groups)
})

test_that("the identifier choice: reference identifiers or the accessions as imported (unchanged)", {
    table <- data.frame(Protein.IDs = "WP_1", gene_name = "a", stringsAsFactors = FALSE)
    object <- methods::new("ProteinQuantitativeData", protein_quant_table = data.frame(Protein.IDs = "WP_1", S1 = 1),
        design_matrix = data.frame(Run = "S1", group = "A"), sample_id = "Run", group_id = "group",
        protein_id_column = "Protein.IDs", protein_id_table = data.frame(Protein.IDs = "WP_1"))
    method <- list(method = "clusterprofiler")
    imported <- resolveProtEnrichAnalysisInputColumns(method, k167Da(table), object)
    expect_identical(imported$idColumn, "Protein.IDs")
    expect_identical(resolveProtEnrichAnalysisInputColumns(method, k167Da(table), object, id_choice = "imported"),
        imported)
    chosen <- resolveProtEnrichAnalysisInputColumns(list(method = "gprofiler2"), k167Da(table), object,
        id_choice = "reference")
    expect_identical(chosen$idColumn, "reference_accession_used")
    expect_identical(chosen$source, "reference_mapping")
})

test_that("the proteomics configuration offers the mapped column only with the gate on", {
    withr::local_options(MultiScholaR.reference_mapping = NULL)
    expect_null(omicGuiConfig("proteomics")$enrichment$id_column$mapped)
    withr::local_options(MultiScholaR.reference_mapping = TRUE)
    expect_identical(omicGuiConfig("proteomics")$enrichment$id_column$mapped, "reference_accession_used")
})

test_that("the shipped GO table holds reference accessions, finds a term for a mapped protein, without network", {
    k167NoNetwork()
    record <- data.frame(field = c("reference_supplied", "uniprot_upid", "mapping_sha256"),
        value = c("FALSE", "UP000000625", strrep("d", 64)), stringsAsFactors = FALSE)
    annotation <- referenceAnnotationTable(data.frame(reference_accession_used = "P0AE06"), record)
    expect_identical(annotation$source, "shipped")
    expect_identical(annotation$go_table_sha256, readReferenceOrganisms()$go_sha256[1])
    expect_true("P0AE06" %in% annotation$table$Entry)
    sets <- buildUniprotGoGeneSets(annotation$table, "none")
    expect_true("P0AE06" %in% sets$term2gene$Entry)
    lost <- attr(annotation$table, "go_ids_not_in_go_db")
    expect_false(is.null(lost$go_db_version))
    expect_identical(lost$total, lost$obsolete + lost$secondary + lost$unknown)
})

test_that("a user-supplied reference is fetched once, cached under the mapping's SHA-256, labelled", {
    raw <- utils::read.delim(referenceOrganismFiles("UP000000625")$go, sep = "\t", quote = "", stringsAsFactors = FALSE)
    raw <- raw[raw$Entry %in% c("P0AE06", "P0AE01"), ]
    raw$From <- raw$Entry
    calls <- 0L
    labels <- character()
    local_mocked_bindings(directUniprotDownload = function(input_tbl, output_path, taxon_id, ..., label = NULL) {
        calls <<- calls + 1L
        labels <<- c(labels, label)
        raw
    })
    project <- tempfile("k167-project-")
    record <- data.frame(field = c("reference_supplied", "uniprot_upid", "uniprot_taxon", "reference", "mapping_sha256"),
        value = c("TRUE", "supplied", "83333", "My K-12", strrep("e", 64)), stringsAsFactors = FALSE)
    mapping <- data.frame(reference_accession_used = c("P0AE06", "P0AE01", NA))
    first <- suppressMessages(referenceAnnotationTable(mapping, record, project_dir = project))
    second <- suppressMessages(referenceAnnotationTable(mapping, record, project_dir = project))
    expect_identical(calls, 1L)
    expect_identical(labels, "reference annotation: My K-12 (supplied)")
    expect_true(file.exists(file.path(project, "reference_mapping", "annotation", strrep("e", 64),
        "uniprot_annotations.RDS")))
    expect_false(is.null(attr(second$table, "go_ids_not_in_go_db")))
    expect_identical(first$go_table_sha256, second$go_table_sha256)
    expect_identical(first$source, "fetched")
})

test_that("the label is appended to the record's query summary; NULL keeps today's text", {
    expect_identical(uniprotLabelSummary("uniprotkb search, organism_id 83333", NULL),
        "uniprotkb search, organism_id 83333")
    expect_identical(uniprotLabelSummary("uniprotkb search, organism_id 83333", "reference annotation: X (UP1)"),
        "uniprotkb search, organism_id 83333; reference annotation: X (UP1)")
})

test_that("the resolver takes the reference branch first and leaves the project's annotation slot unchanged", {
    workflow <- list2env(list(uniprot_dat_cln = data.frame(Entry = "PROJECT")))
    reference <- data.frame(Entry = "P0AE06", gene_names = "acrA")
    resolved <- resolveProtEnrichUniprotAnnotations(workflow, list(source_dir = tempdir(), results_dir = tempdir()),
        referenceAnnotation = reference)
    expect_identical(resolved$source, "reference_mapping")
    expect_identical(resolved$uniprotDatCln$Entry, "P0AE06")
    expect_identical(workflow$uniprot_dat_cln$Entry, "PROJECT")
    without <- resolveProtEnrichUniprotAnnotations(workflow, list(source_dir = tempdir(), results_dir = tempdir()))
    expect_identical(without$source, "workflow_data")
})

test_that("the reference names are declared only when the reference identifiers are used", {
    declared <- diaEnrichmentOperation()$parameters
    plain <- diaEnrichmentDeclaredNames(declared, list())
    expect_false(any(c("reference_mapping_sha256", "id_source", "go_table_sha256") %in% plain))
    used <- diaEnrichmentDeclaredNames(declared, list(id_source = "reference", reference_mapping_sha256 = "x",
        go_table_sha256 = "y"))
    expect_true(all(c("reference_mapping_sha256", "id_source", "go_table_sha256") %in% used))
})

test_that("the integration features gain the mapping columns only with a mapping, by the protein-group rule", {
    ids <- c("ref|WP_1|;ref|WP_2|", "WP_1;WP_4", "P69905")
    annotation <- data.frame(Protein.IDs = ids, gene_name = c("a", "b", "c"), stringsAsFactors = FALSE)
    without <- protIntegrationInputFeatures(ids, annotation)
    with <- protIntegrationInputFeatures(ids, annotation, mapping = k167Mapping)
    expect_identical(with[names(without)], without)
    expect_identical(with$reference_accession_used, c("P1", NA, NA))
    expect_identical(with$map_status, c("mapped", "ambiguous_group", "unmapped"))
    expect_identical(setdiff(names(with), names(without)),
        c("reference_accession_used", "reference_gene", "locus_tag", "map_status"))
})

# --- The run on a strain project (A17): the T7 fixture's proteins (KIT-173) renamed as strain accessions and mapped
# one-to-one to their K-12 accessions; the project's own UniProt annotation is empty by design.
k167EnrichProject <- function(root, id_source, enable_organism_filter = FALSE) {
    withr::local_options(MultiScholaR.reference_mapping = TRUE)
    roles <- goPropReadTsv("t7-fixture-entries.tsv")
    entries <- sort(roles$entry)
    strain <- sprintf("WP_%06d.1", seq_along(entries))
    to_strain <- stats::setNames(strain, entries)
    paths <- goPropPaths(root)
    for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "k167-enrich", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA")
    design <- data.frame(Run = paste0("S", 1:6), group = rep(c("A", "B"), each = 3L), stringsAsFactors = FALSE)
    quant <- goPropQuant(strain, to_strain[roles$entry[roles$role == "p_term"]],
        to_strain[roles$entry[roles$role == "d_term"]])
    object <- methods::new("ProteinQuantitativeData", protein_quant_table = quant, design_matrix = design,
        sample_id = "Run", group_id = "group", protein_id_column = "uniprot_acc",
        protein_id_table = data.frame(uniprot_acc = strain, stringsAsFactors = FALSE),
        args = list(globalParameters = list(workflow_type = "DIA", use_limpa = FALSE, report_template = "DIA_report"),
            deAnalysisParameters = list(formula_string = "~ 0 + group")))
    setDIAStepProvenance(manager, "design_builder")
    manager$saveState("correlation_filtered", object, list(), "fixture")
    manager$.context$experiment_paths <- paths
    mapping <- data.frame(strain_accession = strain, entity_id = strain, reference_accession = entries,
        reference_accession_used = entries, map_status = "one_to_one", stringsAsFactors = FALSE)
    record <- data.frame(field = c("reference", "reference_supplied", "uniprot_upid", "uniprot_taxon", "mapping_sha256"),
        value = c("Escherichia coli K-12 MG1655", "FALSE", "UP000000625", "83333",
            canonicalReferenceMappingSha256(mapping)), stringsAsFactors = FALSE)
    withDIAStateScope(manager, function() {
        saveWorkflowOwnedValue(mapping, "reference_mapping", manager)
        saveWorkflowOwnedValue(record, "reference_mapping_record", manager)
    })
    contrasts <- data.frame(contrasts = "groupB-groupA", full_format = "B_vs_A=groupB-groupA",
        friendly_names = "B_vs_A", stringsAsFactors = FALSE)
    workflow <- list2env(list(state_manager = manager, contrasts_tbl = contrasts, tab_status = list(),
        qc_params = list(protein_qc = list()), config_list = list()))
    parameters <- list(selected_contrast = "B_vs_A", organism_taxid = "999999", up_cutoff = 0, down_cutoff = 0,
        q_cutoff = 0.05, correction_method = "gSCS", enable_organism_filter = enable_organism_filter,
        id_source = id_source)
    run <- NULL
    error <- NULL
    withDIAStateScope(manager, function() {
        da_data <- new.env(parent = emptyenv())
        da_data$current_s4_object <- getDIAStateReference(manager)
        da_data$contrasts_available <- contrasts$contrasts
        manifest <- runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05,
            treat_lfc_cutoff = 0), da_data, workflow, paths, contrasts = contrasts, publish = FALSE)
        produceDIADaPublication(manager, manifest, createEmptyUniprotTable(), 0.05,
            paths[c("da_output_dir", "publication_graphs_dir")], save_assets = FALSE)
        run <<- tryCatch(suppressMessages(runWorkflowEnrichment(workflow, omicGuiConfig("proteomics")$enrichment,
            manifest, parameters, paths, process_fn = processEnrichments)), error = function(e) {
                error <<- e
                NULL
            })
    })
    record <- NULL
    if (!is.null(run)) {
        path <- file.path(root, "analysis_record.json")
        writeWorkflowAnalysisRecord(manager, path)
        record <- jsonlite::fromJSON(path, simplifyVector = FALSE)
    }
    list(run = run, error = error, workflow = workflow, record = record)
}

test_that("a strain whose project annotation is empty runs with the reference identifiers and refuses as imported", {
    skip_if_not_installed("clusterProfiler")
    project <- k167EnrichProject(tempfile("k167-run-"), "reference")
    expect_null(project$error)
    expect_null(project$workflow$uniprot_dat_cln)
    expect_identical(project$run$parameters$id_source, "reference")
    expect_identical(project$run$parameters$go_table_sha256, readReferenceOrganisms()$go_sha256[1])
    expect_false(is.null(project$run$parameters$reference_mapping_sha256))
    tables <- unlist(lapply(project$run$results@enrichment_data, function(contrast) lapply(contrast, function(result)
        if (is.null(result)) 0L else nrow(result@result))))
    expect_gt(sum(tables), 0L)
    expect_gt(project$workflow$enrichment_reference_mapping[[1L]]$mapped, 0L)
    imported <- k167EnrichProject(tempfile("k167-run-"), "imported")
    expect_null(imported$workflow$enrichment_reference_mapping)
    expect_s3_class(imported$error, "workflow_enrichment_annotation")
})

test_that("with the reference identifiers the organism filter is not applied and not declared", {
    skip_if_not_installed("clusterProfiler")
    project <- k167EnrichProject(tempfile("k167-run-"), "reference", enable_organism_filter = TRUE)
    expect_null(project$error)
    expect_false(isTRUE(project$run$organism_filter$filter_applied))
    expect_false(isTRUE(project$run$organism_filter$enabled))
    expect_false(any(c("organism_filter_enabled", "organism_filter_taxon") %in% names(project$run$parameters)))
})

test_that("the analysis record of a reference-identifier run names the three reference values for its replay", {
    skip_if_not_installed("clusterProfiler")
    project <- k167EnrichProject(tempfile("k167-run-"), "reference")
    step <- Filter(function(step) identical(step$operation, "enrichment"), project$record$steps)[[1L]]
    names <- vapply(step$parameters, function(p) p$name, "")
    expect_true(all(c("id_source", "reference_mapping_sha256", "go_table_sha256") %in% names))
    value <- function(name) decodeWorkflowRecordValue(step$parameters[[match(name, names)]]$value)
    expect_identical(value("id_source"), "reference")
    expect_identical(value("go_table_sha256"), readReferenceOrganisms()$go_sha256[1])
})

k167AsRecorded <- function(values) lapply(values, function(value) decodeWorkflowRecordValue(jsonlite::fromJSON(
    jsonlite::toJSON(encodeWorkflowRecordValue(value), auto_unbox = TRUE, digits = I(17), null = "null", na = "null"),
    simplifyVector = FALSE)))

test_that("the replay passes the recorded identifier choice with the gate on and checks the reference values", {
    captured <- NULL
    recorded <- k167AsRecorded(list(taxon_id = "999999", up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05, correction_method = "gSCS",
        id_source = "reference", reference_mapping_sha256 = strrep("a", 64), go_table_sha256 = strrep("b", 64)))
    produced <- recorded
    local_mocked_bindings(produceDIADaPublication = function(...) invisible(NULL),
        runWorkflowEnrichment = function(workflow, config, manifest, parameters, paths, ...) {
            captured <<- list(parameters = parameters, gate = referenceMappingEnabled())
            list(parameters = produced, service = character(), organism_filter = list())
        })
    withr::local_options(MultiScholaR.reference_mapping = NULL)
    contrasts <- data.frame(friendly_names = "B_vs_A", stringsAsFactors = FALSE)
    replayWorkflowEnrichmentStep(NULL, list(state = "enrichment_completed", outcomes = list()), recorded,
        list(parameters = list(qvalue = 0.05)), contrasts, list(source_dir = NULL), tempfile("k167-replay-"))
    expect_identical(captured$parameters$id_source, "reference")
    expect_true(captured$gate)
    expect_false(referenceMappingEnabled())
    produced$go_table_sha256 <- strrep("c", 64)
    expect_error(replayWorkflowEnrichmentStep(NULL, list(state = "enrichment_completed", outcomes = list()), recorded,
        list(parameters = list(qvalue = 0.05)), contrasts, list(source_dir = NULL), tempfile("k167-replay-")),
        class = "workflow_replay_invalid")
})

test_that("a run without the reference identifiers replays without the choice and the gate untouched", {
    captured <- NULL
    recorded <- k167AsRecorded(list(taxon_id = "999999", up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05,
        correction_method = "gSCS"))
    local_mocked_bindings(produceDIADaPublication = function(...) invisible(NULL),
        runWorkflowEnrichment = function(workflow, config, manifest, parameters, paths, ...) {
            captured <<- list(parameters = parameters, gate = referenceMappingEnabled())
            list(parameters = recorded, service = character(), organism_filter = list())
        })
    withr::local_options(MultiScholaR.reference_mapping = NULL)
    replayWorkflowEnrichmentStep(NULL, list(state = "enrichment_completed", outcomes = list()), recorded,
        list(parameters = list(qvalue = 0.05)), data.frame(friendly_names = "B_vs_A"), list(source_dir = NULL),
        tempfile("k167-replay-"))
    expect_false("id_source" %in% names(captured$parameters))
    expect_false(captured$gate)
})

test_that("the 2 x 2 counts use the lists' strict rule and count only tested rows (audit K167-I09, KIT-205)", {
    table <- data.frame(Protein.IDs = c("WP_1", "WP_3", "WP_4", "WP_5"),
        fdr_qvalue = c(0.01, 0.01, NA, 0.01), raw_pvalue = c(0.001, 0.001, NA, 0.001), log2FC = c(1, 1.5, 3, -1),
        stringsAsFactors = FALSE)
    attached <- attachReferenceMapping(k167Da(table), k167Mapping, "Protein.IDs", q_cutoff = 0.05, up_cutoff = 1,
        down_cutoff = 1)
    counts <- attached$counts$B_vs_A
    # WP_1 at log2FC exactly 1 and WP_5 at exactly -1 are not in the lists (`>` and `<`); WP_4 has no p-value.
    expect_identical(counts$untested, 1L)
    expect_identical(counts$groups, 3L)
    expect_identical(c(counts$changed_mapped, counts$changed_unmapped, counts$unchanged_mapped,
        counts$unchanged_unmapped), c(1L, 0L, 1L, 1L))
    expect_false("WP_4" %in% attached$da_results@da_data$B_vs_A$Protein.IDs)
    shown <- referenceMappingContrastTable(counts)
    expect_identical(shown$features, c("mapped", "unmapped", "untested (no p-value in this comparison)"))
    expect_identical(shown$unchanged[3], 1L)
})
