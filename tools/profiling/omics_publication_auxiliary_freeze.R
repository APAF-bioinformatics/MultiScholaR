auxPublicationFreezePaths <- function() {
    list(
        parameters = "tests/testdata/omics-performance/auxiliary/parameters-v1.json",
        sources = "tests/testdata/omics-performance/auxiliary/sources-v1.json",
        responses = "tests/testdata/omics-performance/auxiliary/responses-v1.json",
        splits = "tests/testdata/omics-performance/auxiliary/splits-v1.json",
        surfaces = "tests/testdata/omics-performance/auxiliary/surfaces-v1.json",
        exclusions = "tests/testdata/omics-performance/auxiliary/exclusions-v1.json",
        negatives = paste0(
            "tests/testdata/omics-performance/auxiliary/negative/",
            "contracts-v1.json"
        ),
        workloads = "tests/testdata/omics-performance/auxiliary/workloads",
        truth = "tests/testdata/omics-performance/auxiliary/truth"
    )
}

auxPublicationValidateFreezeAuthorities <- function() {
    paths <- auxPublicationFreezePaths()
    auxPublicationValidateParameters(publicationReadJson(paths$parameters))
    auxPublicationValidateSources(publicationReadJson(paths$sources))
    auxPublicationValidateResponses(publicationReadJson(paths$responses))
    auxPublicationValidateSplits(publicationReadJson(paths$splits))
    auxPublicationValidateSurfaceAuthority(publicationReadJson(paths$surfaces))
    auxPublicationValidateExclusions(publicationReadJson(paths$exclusions))
    auxPublicationValidateNegatives(publicationReadJson(paths$negatives))
    invisible(paths)
}

auxPublicationFreezeOne <- function(
    route_id,
    workload_class,
    output_root,
    contract_root,
    truth_root
) {
    paths <- auxPublicationFreezePaths()
    splits <- publicationReadJson(paths$splits)
    assignment <- auxPublicationFindAssignment(
        splits,
        route_id,
        workload_class
    )
    provisional <- auxPublicationBuildContract(
        route_id,
        workload_class,
        assignment,
        assignment$seed
    )
    auxPublicationValidateContract(provisional)
    workload_id <- provisional$workload_id
    prepared_root <- file.path(output_root, workload_id)
    if (file.exists(prepared_root) || dir.exists(prepared_root)) {
        auxPublicationAbort("auxiliary prepared root already exists")
    }
    dir.create(prepared_root, recursive = TRUE, showWarnings = FALSE)
    payload <- auxPublicationBuildPayload(provisional)
    payload_path <- file.path(prepared_root, provisional$generator$payload_filename)
    auxPublicationWritePayload(payload, payload_path)
    payload_sha256 <- publicationFileDigest(payload_path)
    truth <- auxPublicationBuildTruth(provisional, payload, payload_sha256)
    truth_path <- file.path(
        truth_root,
        paste0(workload_id, ".json")
    )
    if (file.exists(truth_path)) {
        auxPublicationAbort("auxiliary truth already exists")
    }
    publicationWriteJson(truth, truth_path)
    truth_sha256 <- publicationFileDigest(truth_path)
    contract <- auxPublicationBuildContract(
        route_id,
        workload_class,
        assignment,
        assignment$seed,
        payload_sha256,
        truth_sha256
    )
    contract_path <- file.path(
        contract_root,
        paste0(workload_id, ".json")
    )
    if (file.exists(contract_path)) {
        auxPublicationAbort("auxiliary contract already exists")
    }
    publicationWriteJson(contract, contract_path)
    contract <- auxPublicationReadContract(contract_path)
    auxPublicationValidateTruth(truth, contract, payload_path)
    receipt <- list(
        schema = "multischolar.omics_publication_auxiliary_preparation",
        schema_version = .AUX_PUBLICATION_VERSION,
        owner_ticket_id = .AUX_PUBLICATION_OWNER,
        status = "passed",
        workload_id = workload_id,
        contract = list(
            path = contract_path,
            sha256 = publicationFileDigest(contract_path)
        ),
        payload = list(
            path = payload_path,
            sha256 = payload_sha256,
            size_bytes = as.numeric(file.info(payload_path)$size)
        ),
        truth = list(path = truth_path, sha256 = truth_sha256),
        candidate_loaded = FALSE,
        generated_outside_measurement = TRUE,
        publication_authority = FALSE
    )
    receipt_path <- file.path(prepared_root, "preparation-receipt.json")
    publicationWriteJson(receipt, receipt_path)
    list(
        workload_id = workload_id,
        route_id = route_id,
        workload_class = workload_class,
        contract_path = contract_path,
        contract_sha256 = publicationFileDigest(contract_path),
        payload_path = payload_path,
        payload_sha256 = payload_sha256,
        payload_size_bytes = as.numeric(file.info(payload_path)$size),
        truth_path = truth_path,
        truth_sha256 = truth_sha256,
        receipt_path = receipt_path,
        receipt_sha256 = publicationFileDigest(receipt_path),
        candidate_loaded = FALSE,
        publication_authority = FALSE
    )
}

auxPublicationFreeze <- function(
    output_root,
    contract_root = auxPublicationFreezePaths()$workloads,
    truth_root = auxPublicationFreezePaths()$truth,
    workload_classes = c("fixture_correctness", "representative")
) {
    auxPublicationValidateFreezeAuthorities()
    if (file.exists(output_root) || dir.exists(output_root)) {
        auxPublicationAbort("auxiliary freeze output root already exists")
    }
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
    dir.create(contract_root, recursive = TRUE, showWarnings = FALSE)
    dir.create(truth_root, recursive = TRUE, showWarnings = FALSE)
    records <- list()
    for (workload_class in workload_classes) {
        for (route_id in names(auxPublicationRoutes())) {
            records[[length(records) + 1L]] <- auxPublicationFreezeOne(
                route_id,
                workload_class,
                output_root,
                contract_root,
                truth_root
            )
        }
    }
    list(
        schema = "multischolar.omics_publication_auxiliary_freeze",
        schema_version = .AUX_PUBLICATION_VERSION,
        owner_ticket_id = .AUX_PUBLICATION_OWNER,
        status = "passed",
        workload_classes = as.list(workload_classes),
        records = records,
        candidate_loaded = FALSE,
        publication_authority = FALSE
    )
}

auxPublicationFreezeDigestMap <- function(record) {
    stats::setNames(lapply(record$records, function(value) {
        list(
            contract_sha256 = value$contract_sha256,
            payload_sha256 = value$payload_sha256,
            truth_sha256 = value$truth_sha256
        )
    }), vapply(record$records, `[[`, character(1), "workload_id"))
}

auxPublicationCompareFreezes <- function(first, second) {
    first_map <- auxPublicationFreezeDigestMap(first)
    second_map <- auxPublicationFreezeDigestMap(second)
    valid <- identical(names(first_map), names(second_map)) &&
        identical(first_map, second_map)
    if (!valid) auxPublicationAbort("auxiliary freeze determinism differs")
    list(
        schema = "multischolar.omics_publication_auxiliary_determinism",
        schema_version = .AUX_PUBLICATION_VERSION,
        owner_ticket_id = .AUX_PUBLICATION_OWNER,
        status = "passed",
        workload_count = as.integer(length(first_map)),
        digest_map = first_map,
        exact_contract_match = TRUE,
        exact_payload_match = TRUE,
        exact_truth_match = TRUE,
        candidate_loaded = FALSE,
        publication_authority = FALSE
    )
}
