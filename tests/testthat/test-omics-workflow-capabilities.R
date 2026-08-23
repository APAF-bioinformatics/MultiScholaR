.omics_capability_path <- function() {
  testthat::test_path("..", "testdata", "omics-capabilities.json")
}

.read_omics_capabilities <- function() {
  jsonlite::read_json(.omics_capability_path(), simplifyVector = FALSE)
}

.omics_repo_path <- function(path) {
  file.path(testthat::test_path("..", ".."), path)
}

.flatten_omics_capabilities <- function(inventory) {
  do.call(c, lapply(inventory$formats, `[[`, "capabilities"))
}

.capability_key <- function(identity, fields) {
  paste(vapply(fields, function(field) identity[[field]], character(1)), collapse = "\x1f")
}

.resolve_inventory_backend <- function(inventory, identity, requested_backend) {
  stopifnot(requested_backend %in% inventory$backend_contract$allowed_requests)

  if (identical(requested_backend, "memory")) {
    return(list(effective_backend = "memory", reason_code = "explicit_memory"))
  }

  fields <- unlist(inventory$capability_key_fields, use.names = FALSE)
  requested_key <- .capability_key(identity, fields)
  capabilities <- .flatten_omics_capabilities(inventory)
  matches <- vapply(
    capabilities,
    function(capability) identical(.capability_key(capability$identity, fields), requested_key),
    logical(1)
  )

  if (sum(matches) != 1L) {
    policy <- inventory$backend_contract$unknown_tuple_policy[[requested_backend]]
    if (startsWith(policy, "error_")) {
      return(list(effective_backend = NULL, reason_code = policy))
    }
    return(list(effective_backend = policy, reason_code = "unknown_capability"))
  }

  capability <- capabilities[[which(matches)]]
  if (identical(requested_backend, "artifact") && !isTRUE(capability$artifact_eligible)) {
    return(list(effective_backend = NULL, reason_code = "error_not_certified"))
  }
  if (identical(requested_backend, "auto") && !isTRUE(capability$auto_eligible)) {
    return(list(effective_backend = "memory", reason_code = "not_promoted"))
  }

  list(
    effective_backend = "artifact",
    reason_code = if (identical(requested_backend, "auto")) "auto_promoted" else "artifact_certified"
  )
}

test_that("all-omics capability inventory has a stable additive schema", {
  inventory <- .read_omics_capabilities()
  formats <- inventory$formats
  capabilities <- .flatten_omics_capabilities(inventory)

  expect_match(inventory$schema_version, "^[0-9]+\\.[0-9]+\\.[0-9]+$")
  expect_match(inventory$inventory_version, "^[0-9]{4}-[0-9]{2}-[0-9]{2}\\.[0-9]+$")
  expect_identical(
    unlist(inventory$session_identity_fields, use.names = FALSE),
    c(
      "project_id", "omic_type", "omic_label", "workflow_id",
      "workflow_type", "workflow_slug", "input_format", "data_level",
      "acquisition_mode"
    )
  )
  expect_identical(
    unlist(inventory$capability_key_fields, use.names = FALSE),
    c(
      "omic_type", "workflow_id", "workflow_type", "workflow_slug",
      "input_format", "data_level", "acquisition_mode"
    )
  )

  format_ids <- vapply(formats, `[[`, character(1), "format_id")
  capability_ids <- vapply(capabilities, `[[`, character(1), "capability_id")
  expect_identical(anyDuplicated(format_ids), 0L)
  expect_identical(anyDuplicated(capability_ids), 0L)
  expect_length(formats, 16L)
  expect_length(capabilities, 11L)

  allowed_statuses <- unlist(inventory$allowed_support_statuses, use.names = FALSE)
  statuses <- vapply(formats, `[[`, character(1), "support_status")
  expect_true(all(statuses %in% allowed_statuses))
  expect_setequal(
    allowed_statuses,
    c(
      "scientifically_supported", "reader_characterized", "detection_only",
      "advertised_unverified", "unsupported"
    )
  )

  key_fields <- unlist(inventory$capability_key_fields, use.names = FALSE)
  capability_keys <- vapply(capabilities, function(capability) {
    expect_setequal(names(capability$identity), key_fields)
    expect_match(capability$capability_version, "^[0-9]+\\.[0-9]+\\.[0-9]+$")
    expect_true(all(vapply(capability$identity, function(value) {
      is.character(value) && length(value) == 1L && nzchar(value)
    }, logical(1))))
    .capability_key(capability$identity, key_fields)
  }, character(1))
  expect_identical(anyDuplicated(capability_keys), 0L)
})

test_that("inventory covers every source-owned import UI choice exactly", {
  inventory <- .read_omics_capabilities()
  format_ids <- vapply(inventory$formats, `[[`, character(1), "format_id")
  referenced_format_ids <- character()

  for (surface in inventory$ui_surfaces) {
    choice_fn <- get(surface$choice_function, envir = asNamespace("MultiScholaR"))
    actual <- choice_fn()
    expected_values <- vapply(surface$choices, `[[`, character(1), "value")
    expected_labels <- vapply(surface$choices, `[[`, character(1), "label")

    expect_identical(unname(actual), expected_values, info = surface$omic_type)
    expect_identical(names(actual), expected_labels, info = surface$omic_type)

    surface_format_ids <- vapply(surface$choices, function(choice) {
      choice$inventory_format_id %||% NA_character_
    }, character(1))
    expect_true(all(is.na(surface_format_ids) | surface_format_ids %in% format_ids))
    referenced_format_ids <- c(referenced_format_ids, stats::na.omit(surface_format_ids))
  }

  expect_setequal(referenced_format_ids, format_ids)
  expect_identical(anyDuplicated(referenced_format_ids), 0L)
  expect_setequal(
    vapply(inventory$ui_surfaces, `[[`, character(1), "omic_type"),
    c("proteomics", "metabolomics", "lipidomics")
  )
})

test_that("parser ownership and evidence tiers do not overstate support", {
  inventory <- .read_omics_capabilities()
  namespace <- asNamespace("MultiScholaR")

  for (format in inventory$formats) {
    expect_true(exists(format$detector_owner, envir = namespace, mode = "function"))
    if (!is.null(format$parser_owner)) {
      expect_true(exists(format$parser_owner, envir = namespace, mode = "function"))
    }
    if (!is.null(format$fallback_parser_owner)) {
      expect_true(exists(format$fallback_parser_owner, envir = namespace, mode = "function"))
    }

    required_tests <- unlist(format$required_tests, use.names = FALSE)
    fixture_refs <- unlist(format$fixture_refs, use.names = FALSE)
    expect_true(length(required_tests) > 0L, info = format$format_id)
    expect_true(all(file.exists(.omics_repo_path(required_tests))), info = format$format_id)
    expect_true(all(file.exists(.omics_repo_path(fixture_refs))), info = format$format_id)

    if (identical(format$support_status, "scientifically_supported")) {
      expect_true(startsWith(format$fixture_requirement, "direct_file"))
      expect_gt(length(fixture_refs), 0L)
      expect_true(all(startsWith(fixture_refs, "tests/testdata/")))
      expect_gt(length(format$capabilities), 0L)
      expect_true(all(vapply(format$capabilities, `[[`, logical(1), "complete_workflow")))
      expect_true(all(vapply(format$capabilities, function(capability) {
        length(capability$e2e_lanes) > 0L
      }, logical(1))))
    }

    if (identical(format$support_status, "reader_characterized")) {
      expect_false(is.null(format$parser_owner))
      expect_true(all(vapply(format$capabilities, function(capability) {
        !isTRUE(capability$complete_workflow) && length(capability$e2e_lanes) == 0L
      }, logical(1))))
    }

    if (identical(format$support_status, "detection_only")) {
      expect_null(format$parser_owner)
      expect_null(format$fallback_parser_owner)
      expect_length(format$capabilities, 0L)
    }
  }
})

test_that("capabilities preserve current analytical routing and coverage", {
  inventory <- .read_omics_capabilities()
  capabilities <- .flatten_omics_capabilities(inventory)
  e2e_manifest <- jsonlite::read_json(
    testthat::test_path("..", "testdata", "e2e", "manifest.json"),
    simplifyVector = FALSE
  )
  e2e_lane_ids <- vapply(e2e_manifest$lanes, `[[`, character(1), "lane_id")

  for (capability in capabilities) {
    expect_true(all(unlist(capability$module_ci_lanes, use.names = FALSE) %in%
      c("proteomics", "metabolomics", "lipidomics")))
    expect_true(all(unlist(capability$e2e_lanes, use.names = FALSE) %in% e2e_lane_ids))
    expect_false(isTRUE(capability$auto_eligible))
    if (identical(capability$capability_id, "proteomics.diann.peptide.dia.v1")) {
      expect_true(isTRUE(capability$artifact_eligible))
      expect_identical(capability$maximum_artifact_rollout, "dual_write")
    } else {
      expect_false(isTRUE(capability$artifact_eligible))
      expect_null(capability$maximum_artifact_rollout)
    }
  }

  spectronaut <- Filter(
    function(capability) identical(capability$identity$input_format, "spectronaut"),
    capabilities
  )
  expect_length(spectronaut, 2L)
  expect_setequal(vapply(spectronaut, function(x) x$identity$data_level, character(1)), c("protein", "peptide"))
  expect_true(all(vapply(spectronaut, function(x) identical(x$identity$workflow_type, "LFQ"), logical(1))))
  expect_true(all(vapply(spectronaut, function(x) identical(x$identity$acquisition_mode, "dia"), logical(1))))

  assay_omics <- Filter(
    function(capability) capability$identity$omic_type %in% c("metabolomics", "lipidomics"),
    capabilities
  )
  expect_true(all(vapply(assay_omics, function(x) {
    identical(x$identity$acquisition_mode, "not_recorded")
  }, logical(1))))
})

test_that("exact capability matching fails closed for artifact requests", {
  inventory <- .read_omics_capabilities()
  known <- .flatten_omics_capabilities(inventory)[[1]]$identity

  expect_identical(
    .resolve_inventory_backend(inventory, known, "memory"),
    list(effective_backend = "memory", reason_code = "explicit_memory")
  )
  expect_identical(
    .resolve_inventory_backend(inventory, known, "artifact"),
    list(effective_backend = "artifact", reason_code = "artifact_certified")
  )
  expect_identical(
    .resolve_inventory_backend(inventory, known, "auto"),
    list(effective_backend = "memory", reason_code = "not_promoted")
  )

  unknown <- known
  unknown$acquisition_mode <- "unregistered_mode"
  expect_identical(
    .resolve_inventory_backend(inventory, unknown, "artifact"),
    list(effective_backend = NULL, reason_code = "error_unknown_capability")
  )
  expect_identical(
    .resolve_inventory_backend(inventory, unknown, "auto"),
    list(effective_backend = "memory", reason_code = "unknown_capability")
  )
})

test_that("placeholder and auxiliary omics surfaces are not promoted as workflows", {
  inventory <- .read_omics_capabilities()
  surfaces <- inventory$non_workflow_surfaces

  expect_setequal(
    vapply(surfaces, `[[`, character(1), "surface_id"),
    c("transcriptomics.ui", "integration.ui", "phosphosite.api", "multiomics.api")
  )
  expect_setequal(
    vapply(surfaces, `[[`, character(1), "surface_type"),
    c("ui_placeholder", "auxiliary_exported_api")
  )
  for (surface in surfaces) {
    expect_true(
      all(file.exists(.omics_repo_path(unlist(surface$owners, use.names = FALSE)))),
      info = surface$surface_id
    )
    expect_true(surface$programme_treatment %in% c(
      "excluded_until_implemented",
      "compatibility_and_performance_only"
    ))
  }
})
