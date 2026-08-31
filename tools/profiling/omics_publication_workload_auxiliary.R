.AUX_PUBLICATION_ADAPTER_ID <- "multischolar.publication.auxiliary"
.AUX_PUBLICATION_ADAPTER_VERSION <- "1.0.0"

auxPublicationReadContract <- function(path) {
    contract <- publicationReadJson(path)
    auxPublicationValidateContract(contract)
    contract
}

auxPublicationBuildPayload <- function(contract) {
    auxPublicationValidateContract(contract)
    parameters <- publicationReadJson(contract$parameter_authority$path)
    responses <- publicationReadJson(contract$response_authority$path)
    auxPublicationValidateParameters(parameters)
    auxPublicationValidateResponses(responses)
    dimensions <- contract$dimensions
    route_id <- contract$route$route_id
    seed <- if (is.null(contract$rng$seed)) 680001L else {
        as.integer(contract$rng$seed)
    }
    if (identical(route_id, "phosphosite_stages")) {
        plan <- auxPublicationPhosphositeFeaturePlan(
            parameters,
            as.integer(dimensions$feature_count),
            as.integer(dimensions$measured_sample_count),
            seed
        )
        payload <- c(list(
            schema = "multischolar.omics_publication_auxiliary_payload",
            schema_version = .AUX_PUBLICATION_VERSION,
            route_id = route_id,
            workload_id = contract$workload_id,
            dimensions = dimensions,
            responses = list()
        ), plan)
        auxPublicationValidatePayload(payload, contract)
        return(payload)
    }
    plan <- auxPublicationMultiomicsPlan(
        parameters,
        as.integer(dimensions$feature_count),
        as.integer(dimensions$sample_count),
        seed
    )
    if (identical(route_id, "mofa_weights")) {
        plan$assay_data <- list()
        plan$mapping_table <- plan$mapping_table[0L, , drop = FALSE]
        plan$measured_sample_ids <- character()
    } else if (identical(route_id, "stringdb_rank")) {
        plan$model_double$weights <- list()
        plan$assay_data <- list()
        plan$mapping_table <- plan$mapping_table[0L, , drop = FALSE]
        plan$measured_sample_ids <- character()
    } else {
        plan$model_double$weights <- list()
    }
    plan$responses <- auxPublicationResponsePlan(
        parameters,
        responses,
        unlist(contract$route$service_ids, use.names = FALSE),
        as.integer(dimensions$response_row_count)
    )
    payload <- c(list(
        schema = "multischolar.omics_publication_auxiliary_payload",
        schema_version = .AUX_PUBLICATION_VERSION,
        route_id = route_id,
        workload_id = contract$workload_id,
        dimensions = dimensions
    ), plan)
    auxPublicationValidatePayload(payload, contract)
    payload
}

auxPublicationValidatePhosphositePayload <- function(payload, contract) {
    dimensions <- contract$dimensions
    required_evidence <- c(
        "leading_proteins", "phospho_sty_probabilities", "phospho_sty",
        "sequence", "experiment"
    )
    sample_columns <- grep(
        "^reporter intensity corrected_[0-9]+$",
        names(payload$evidence),
        value = TRUE
    )
    categories <- c("single", "ambiguous", "missing", "contaminant", "reverse")
    valid <- is.data.frame(payload$evidence) &&
        is.data.frame(payload$proteins) && is.data.frame(payload$row_plan) &&
        all(required_evidence %in% names(payload$evidence)) &&
        identical(nrow(payload$evidence), as.integer(dimensions$feature_count)) &&
        identical(nrow(payload$row_plan), as.integer(dimensions$feature_count)) &&
        identical(
            length(sample_columns),
            as.integer(dimensions$measured_sample_count)
        ) &&
        identical(payload$row_plan$source_index, seq_len(nrow(payload$row_plan))) &&
        all(payload$row_plan$category_id %in% categories) &&
        all(payload$evidence$phospho_sty %in% c(1L, 2L)) &&
        all(grepl(
            "^[A-Z]+\\([0-9]+[.][0-9]+\\)[A-Z]",
            payload$evidence$phospho_sty_probabilities
        )) &&
        all(grepl("^[A-Z]+$", payload$evidence$sequence)) &&
        !anyNA(payload$proteins$seq) &&
        all(grepl("^[A-Z]+$", payload$proteins$seq)) &&
        !anyDuplicated(payload$proteins$uniprot_acc)
    if (!valid) auxPublicationAbort("auxiliary phosphosite payload differs")
    invisible(payload)
}

auxPublicationValidateMultiomicsPayload <- function(payload, contract) {
    dimensions <- contract$dimensions
    expected_layers <- c("proteome", "metabolome_lc", "metabolome_gc")
    layer_counts <- table(factor(
        payload$feature_plan$layer,
        levels = expected_layers
    ))
    sample_counts <- vapply(payload$sample_registry, length, integer(1))
    assay_rows <- sum(vapply(payload$assay_data, nrow, integer(1)))
    expected_assay_rows <- sum(payload$weights$view %in% c(
        "metabolome_lc",
        "metabolome_gc"
    ))
    route_id <- contract$route$route_id
    route_payload_valid <- if (identical(route_id, "mofa_weights")) {
        identical(names(payload$model_double$weights), expected_layers) &&
            !length(payload$assay_data) && !nrow(payload$mapping_table) &&
            !length(payload$measured_sample_ids)
    } else if (route_id %in% c(
        "metabolite_enrichment",
        "metabolite_pathway"
    )) {
        !length(payload$model_double$weights) &&
            identical(
                names(payload$assay_data),
                c("metabolome_lc", "metabolome_gc")
            ) && identical(length(payload$measured_sample_ids), 1L)
    } else {
        !length(payload$model_double$weights) &&
            !length(payload$assay_data) && !nrow(payload$mapping_table) &&
            !length(payload$measured_sample_ids)
    }
    valid <- route_payload_valid &&
        identical(names(payload$sample_registry), expected_layers) &&
        identical(sum(layer_counts), as.integer(dimensions$feature_count)) &&
        identical(nrow(payload$weights), as.integer(dimensions$feature_count)) &&
        identical(nrow(payload$feature_plan), as.integer(dimensions$feature_count)) &&
        all(sample_counts == as.integer(dimensions$sample_count)) &&
        identical(assay_rows, as.integer(
            if (route_id %in% c(
                "metabolite_enrichment",
                "metabolite_pathway"
            )) expected_assay_rows else 0L
        )) &&
        identical(
            length(payload$measured_sample_ids),
            as.integer(dimensions$measured_sample_count)
        ) &&
        all(c("feature", "view", "factor", "value") %in% names(payload$weights)) &&
        all(payload$weights$factor == "Factor1") &&
        all(payload$weights$view %in% expected_layers)
    if (!valid) auxPublicationAbort("auxiliary multiomics payload differs")
    parameters <- publicationReadJson(contract$parameter_authority$path)
    responses <- publicationReadJson(contract$response_authority$path)
    multiplicity <- as.integer(auxPublicationParameter(
        parameters,
        "multiomics_response_identifier_multiplicity"
    ))
    pathway_count <- as.integer(auxPublicationParameter(
        parameters,
        "multiomics_response_pathway_count"
    ))
    expected_services <- unlist(contract$route$service_ids, use.names = FALSE)
    if (!identical(names(payload$responses), expected_services)) {
        auxPublicationAbort("auxiliary payload response services differ")
    }
    lapply(expected_services, function(service_id) {
        definition <- auxPublicationResponseDefinition(responses, service_id)
        auxPublicationValidateExpandedResponse(
            payload$responses[[service_id]],
            definition,
            as.integer(dimensions$response_row_count),
            multiplicity,
            pathway_count
        )
    })
    invisible(payload)
}

auxPublicationValidatePayload <- function(payload, contract) {
    if (identical(contract$route$route_id, "phosphosite_stages")) {
        return(auxPublicationValidatePhosphositePayload(payload, contract))
    }
    auxPublicationValidateMultiomicsPayload(payload, contract)
}

auxPublicationWritePayload <- function(payload, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- tempfile("auxiliary-payload-", tmpdir = dirname(path))
    on.exit(unlink(temporary, force = TRUE), add = TRUE)
    saveRDS(payload, temporary, version = 3L, compress = FALSE)
    if (!file.rename(temporary, path)) {
        auxPublicationAbort("could not atomically publish auxiliary payload")
    }
    invisible(path)
}

auxPublicationLoadPayload <- function(path, contract) {
    if (!file.exists(path) || !identical(
        publicationFileDigest(path),
        contract$expected_digests$payload_sha256
    )) {
        auxPublicationAbort("auxiliary payload binding differs")
    }
    payload <- readRDS(path)
    publicationRequireNames(
        payload,
        if (identical(contract$route$route_id, "phosphosite_stages")) {
            c(
                "schema", "schema_version", "route_id", "workload_id",
                "dimensions", "responses", "evidence", "proteins", "row_plan"
            )
        } else {
            c(
                "schema", "schema_version", "route_id", "workload_id",
                "dimensions", "model_double", "weights", "assay_data",
                "mapping_table", "sample_registry", "measured_sample_ids",
                "feature_plan", "responses"
            )
        },
        "Auxiliary payload"
    )
    valid <- identical(
        payload$schema,
        "multischolar.omics_publication_auxiliary_payload"
    ) && identical(payload$schema_version, .AUX_PUBLICATION_VERSION) &&
        identical(payload$route_id, contract$route$route_id) &&
        identical(payload$workload_id, contract$workload_id) &&
        identical(
            publicationObjectDigest(payload$dimensions),
            publicationObjectDigest(contract$dimensions)
        )
    if (!valid) auxPublicationAbort("auxiliary payload differs")
    auxPublicationValidatePayload(payload, contract)
    payload
}

auxPublicationFunctionWithOverrides <- function(fun, replacements) {
    clone <- fun
    environment(clone) <- list2env(replacements, parent = environment(fun))
    clone
}

auxPublicationExport <- function(name) {
    getExportedValue("MultiScholaR", name)
}

auxPublicationCallTidy <- function(fun, args, tidy_args = list()) {
    call <- rlang::call2(fun, !!!args, !!!tidy_args)
    rlang::eval_tidy(call)
}

auxPublicationRunPhosphosite <- function(payload) {
    previous_plan <- future::plan()
    on.exit(future::plan(previous_plan), add = TRUE)
    future::plan(future::sequential)
    clean <- auxPublicationCallTidy(
        auxPublicationExport("addColumnsToEvidenceTbl"),
        list(payload$evidence),
        list(rlang::sym("phospho_sty_probabilities"))
    )
    accession <- auxPublicationCallTidy(
        auxPublicationExport("chooseBestPhosphositeAccession"),
        list(clean, payload$proteins),
        list(
            rlang::sym("leading_proteins"),
            rlang::sym("evidence_id")
        )
    )
    abundance <- auxPublicationExport("removePeptidesWithoutAbundances")(
        clean,
        "reporter intensity corrected"
    )
    probabilities <- auxPublicationCallTidy(
        auxPublicationExport("filterPeptideAndExtractProbabilities"),
        list(
            abundance,
            accession,
            "reporter intensity corrected"
        ),
        list(
            accession_col = rlang::sym("leading_proteins"),
            phospho_site_prob_col = rlang::sym(
                "phospho_sty_probabilities"
            ),
            num_phospho_site_col = rlang::sym("phospho_sty")
        )
    )
    peptide_bounds <- auxPublicationExport("addPeptideStartAndEnd")(
        probabilities,
        payload$proteins
    )
    positions <- auxPublicationExport("addPhosphositesPositionsString")(
        peptide_bounds
    )
    xmers <- auxPublicationExport("addXMerStrings")(positions, 7L)
    filtered <- auxPublicationCallTidy(
        auxPublicationExport("filterByScoreAndGetSimilarPeptides"),
        list(xmers, 0.75, 0.5),
        list(num_phospho_site_col = rlang::sym("phospho_sty"))
    )
    long <- auxPublicationCallTidy(
        auxPublicationExport("allPhosphositesPivotLonger"),
        list(
            filtered,
            "experiment",
            "reporter intensity corrected",
            "_\\d+",
            "_(\\d+)"
        ),
        list(
            phospho_site_prob_col = rlang::sym(
                "phospho_sty_probabilities"
            ),
            num_phospho_site_col = rlang::sym("phospho_sty")
        )
    )
    paralog <- auxPublicationCallTidy(
        auxPublicationExport("groupParalogPeptides"),
        list(long, "experiment"),
        list(
            phospho_site_prob_col = rlang::sym(
                "phospho_sty_probabilities"
            ),
            num_phospho_site_col = rlang::sym("phospho_sty")
        )
    )
    wide <- auxPublicationCallTidy(
        auxPublicationExport("allPhosphositesPivotWider"),
        list(paralog, "experiment"),
        list(
            phospho_site_prob_col = rlang::sym(
                "phospho_sty_probabilities"
            ),
            num_phospho_site_col = rlang::sym("phospho_sty")
        )
    )
    summarised_long <- auxPublicationExport(
        "uniquePhosphositesSummariseLongList"
    )(paralog, "experiment")
    summarised_wide <- auxPublicationExport(
        "uniquePhosphositesSummariseWideList"
    )(summarised_long, "experiment")
    list(
        stages = list(
            cleaned = nrow(clean),
            accession = nrow(accession),
            abundance = nrow(abundance),
            probabilities = nrow(probabilities),
            peptide_bounds = nrow(peptide_bounds),
            positions = nrow(positions),
            xmers = nrow(xmers),
            filtered = nrow(filtered),
            long = nrow(long),
            paralog = nrow(paralog),
            wide = nrow(wide)
        ),
        retained = list(
            summarised_wide_list = summarised_wide,
            summarised_long_list = summarised_long,
            all_phos_sites_wide = wide,
            all_phos_sites_long = long
        )
    )
}

auxPublicationPlotFunction <- function(output_root) {
    fun <- auxPublicationExport("plotMofaWeights")
    body_text <- paste(deparse(body(fun)), collapse = "\n")
    body_text <- gsub(
        "MOFA2::get_weights",
        "auxiliaryGetWeights",
        body_text,
        fixed = TRUE
    )
    body(fun) <- parse(text = body_text)[[1L]]
    auxPublicationFunctionWithOverrides(fun, list(
        auxiliaryGetWeights = function(model) model$weights,
        getProjectPaths = function(...) list(mofa_plots_dir = output_root),
        omic_type = "multiomics",
        experiment_label = "publication_auxiliary"
    ))
}

auxPublicationRunMofaWeights <- function(payload, output_root) {
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
    plot_fun <- auxPublicationPlotFunction(output_root)
    plots <- lapply(names(payload$model_double$weights), function(view) {
        plot_fun(payload$model_double, view, factor_level = Factor1)
    })
    names(plots) <- names(payload$model_double$weights)
    list(
        stages = list(
            view_count = as.integer(length(plots)),
            ranked_feature_count = as.integer(nrow(payload$weights)),
            plotted_feature_count = as.integer(length(plots) * 20L),
            file_count = as.integer(length(list.files(output_root)))
        ),
        retained = plots
    )
}

auxPublicationMetaboliteObject <- function(payload) {
    design <- data.frame(
        sample_id = payload$measured_sample_ids,
        group = "synthetic",
        stringsAsFactors = FALSE
    )
    methods::new(
        "MetaboliteAssayData",
        metabolite_data = payload$assay_data,
        metabolite_id_column = "metabolite",
        annotation_id_column = "database_identifier",
        database_identifier_type = "CHEBI",
        internal_standard_regex = NA_character_,
        design_matrix = design,
        sample_id = "sample_id",
        group_id = "group",
        technical_replicate_id = NA_character_,
        args = list()
    )
}

auxPublicationRunMetaboliteEnrichment <- function(payload, output_root) {
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
    object <- auxPublicationMetaboliteObject(payload)
    observed <- new.env(parent = emptyenv())
    observed$service <- character()
    observed$ranked_list_length <- integer()
    response <- function(service_id) {
        function(ranked_list, ...) {
            observed$service <- c(observed$service, service_id)
            observed$ranked_list_length <- c(
                observed$ranked_list_length,
                length(ranked_list)
            )
            payload$responses[[service_id]]
        }
    }
    fun <- auxPublicationFunctionWithOverrides(
        auxPublicationExport("runMetabolomicsEnrichmentAnalysis"),
        list(
            getProjectPaths = function(...) list(
                integration_enrichment_plots_dir = output_root
            ),
            runKeggEnrichment = response("kegg"),
            runReactomeEnrichment = response("reactome")
        )
    )
    results <- lapply(c("metabolome_lc", "metabolome_gc"), function(view) {
        fun(
            weights = payload$weights,
            metabolomics_obj = object,
            mapping_table = payload$mapping_table,
            project_dirs = list(),
            omic_type = "multiomics",
            experiment_label = "publication_auxiliary",
            assay_name = view,
            kegg_species_code = "synthetic",
            reactome_organism = "Synthetic organism"
        )
    })
    names(results) <- c("metabolome_lc", "metabolome_gc")
    list(
        stages = list(
            input_feature_count = as.integer(nrow(payload$weights)),
            response_service_order = as.list(observed$service),
            ranked_list_lengths = as.list(observed$ranked_list_length),
            result_row_count = as.integer(sum(vapply(results, nrow, integer(1)))),
            file_count = as.integer(length(list.files(output_root)))
        ),
        retained = results
    )
}

auxPublicationRunMetabolitePathway <- function(payload, output_root) {
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
    object <- auxPublicationMetaboliteObject(payload)
    response <- payload$responses$metabolomics_pathway
    calls <- new.env(parent = emptyenv())
    calls$assay <- character()
    fun <- auxPublicationFunctionWithOverrides(
        auxPublicationExport("runMetabolomicsPathwayEnrichment"),
        list(
            getProjectPaths = function(...) list(
                integration_enrichment_plots_dir = output_root
            ),
            runMetabolomicsEnrichmentAnalysis = function(..., assay_name) {
                calls$assay <- c(calls$assay, assay_name)
                response
            }
        )
    )
    result <- fun(
        weights = payload$weights,
        metabolomics_obj = object,
        mapping_table = payload$mapping_table,
        project_dirs = list(),
        omic_type = "multiomics",
        experiment_label = "publication_auxiliary",
        kegg_species_code = "synthetic",
        reactome_organism = "Synthetic organism"
    )
    list(
        stages = list(
            input_feature_count = as.integer(nrow(payload$weights)),
            assay_call_order = as.list(calls$assay),
            result_row_count = as.integer(nrow(result)),
            mapped_name_count = as.integer(sum(
                !is.na(result$mappedNames) & nzchar(result$mappedNames)
            )),
            file_count = as.integer(length(list.files(output_root)))
        ),
        retained = result
    )
}

auxPublicationRunStringDbRank <- function(payload, output_root) {
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
    response <- payload$responses$stringdb
    fun <- auxPublicationFunctionWithOverrides(
        auxPublicationExport("runOneStringDbRankEnrichmentMofa"),
        list(
            submitStringDBEnrichment = function(...) list(
                job_id = "synthetic-local-job",
                api_key = "synthetic-local-key",
                submission_response = list(status = "submitted")
            ),
            retrieveStringDBEnrichmentResults = function(...) list(
                enrichment_data = response,
                page_url = "https://invalid.example/local/page",
                download_url = "https://invalid.example/local/result",
                graph_url = "https://invalid.example/local/graph",
                graph_image_content = as.raw(c(137L, 80L, 78L, 71L))
            )
        )
    )
    input <- payload$weights[, c("feature", "value"), drop = FALSE]
    result <- fun(
        input_table = input,
        identifier_column_name = "feature",
        value_column_name = "value",
        result_label = "synthetic_local",
        results_dir = output_root,
        api_key = "synthetic-local-key",
        polling_interval_seconds = 1L,
        max_polling_attempts = 1L
    )
    list(
        stages = list(
            ranked_identifier_count = as.integer(nrow(input)),
            result_row_count = as.integer(nrow(result)),
            file_count = as.integer(length(list.files(output_root)))
        ),
        retained = result
    )
}

auxPublicationRunPayload <- function(contract, payload, output_root) {
    switch(
        contract$route$route_id,
        phosphosite_stages = auxPublicationRunPhosphosite(payload),
        mofa_weights = auxPublicationRunMofaWeights(payload, output_root),
        metabolite_enrichment = auxPublicationRunMetaboliteEnrichment(
            payload,
            output_root
        ),
        metabolite_pathway = auxPublicationRunMetabolitePathway(
            payload,
            output_root
        ),
        stringdb_rank = auxPublicationRunStringDbRank(payload, output_root),
        auxPublicationAbort("auxiliary workload route is unsupported")
    )
}

auxPublicationResultEvidence <- function(result) {
    list(
        stages = result$stages,
        retained_sha256 = digest::digest(
            result$retained,
            algo = "sha256",
            serialize = TRUE
        ),
        retained_class = as.list(class(result$retained)),
        promotion_authority = FALSE
    )
}

omicsWorkloadAdapter <- function() {
    list(
        adapter_id = .AUX_PUBLICATION_ADAPTER_ID,
        adapter_version = .AUX_PUBLICATION_ADAPTER_VERSION,
        supported_omics = "auxiliary",
        prepare = function(context) {
            payload_path <- file.path(context$run_dir, "prepared", "payload.rds")
            auxPublicationLoadPayload(payload_path, context$contract)
            list(
                payload_root = dirname(payload_path),
                payload_path = payload_path,
                truth_path = file.path(dirname(payload_path), "truth.json")
            )
        },
        run = function(context) {
            payload <- auxPublicationLoadPayload(
                context$payload_path,
                context$contract
            )
            result <- auxPublicationRunPayload(
                context$contract,
                payload,
                file.path(context$run_dir, "output")
            )
            truth <- publicationReadJson(context$truth_path)
            auxPublicationValidateTruth(
                truth,
                context$contract,
                context$payload_path
            )
            auxPublicationValidateRunResult(result, truth, context$contract)
            file_count <- if (is.null(result$stages$file_count)) {
                0L
            } else {
                result$stages$file_count
            }
            list(
                status = "passed",
                workflow_evidence = auxPublicationResultEvidence(result),
                query_evidence = list(),
                session_evidence = list(
                    status = "memory_only",
                    resource_count = 0L
                ),
                report_evidence = list(
                    status = "local_products",
                    file_count = file_count
                ),
                retained = result$retained
            )
        }
    )
}
