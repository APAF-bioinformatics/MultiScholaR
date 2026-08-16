.MODULE_CI_VALID_OMICS <- c("all", "proteomics", "metabolomics", "lipidomics")
.MODULE_CI_DEFAULT_RUNTIME_CLASSES <- c(
  "unit-contract",
  "module-browser",
  "module-artifact",
  "workflow-e2e",
  "release-full"
)

module_ci_coalesce <- function(x, y) {
  if (is.null(x)) y else x
}

module_ci_fixture_root <- function() {
  testthat::test_path("..", "testdata", "mci")
}

module_ci_manifest_path <- function() {
  file.path(module_ci_fixture_root(), "manifest.json")
}

module_ci_required_scenario_fields <- function() {
  c(
    "scenario_id",
    "ticket_id",
    "item_id",
    "omic",
    "module_family",
    "runtime",
    "ci_lane",
    "fixture",
    "parameters",
    "expected_state",
    "expected_artifacts",
    "required_tests"
  )
}

read_module_ci_manifest <- function(
    manifest_path = module_ci_manifest_path(),
    fixture_root = dirname(manifest_path),
    validate = TRUE
) {
  if (!file.exists(manifest_path)) {
    stop(sprintf("module CI manifest not found: %s", manifest_path), call. = FALSE)
  }

  manifest <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
  if (isTRUE(validate)) {
    validate_module_ci_manifest(manifest, fixture_root = fixture_root, manifest_path = manifest_path)
  }

  scenarios <- manifest$scenarios
  names(scenarios) <- vapply(scenarios, `[[`, character(1), "scenario_id")
  manifest$scenarios <- scenarios
  manifest
}

validate_module_ci_manifest <- function(
    manifest,
    fixture_root = module_ci_fixture_root(),
    manifest_path = module_ci_manifest_path()
) {
  errors <- character()

  add_error <- function(message) {
    errors <<- c(errors, message)
  }

  if (is.null(manifest$schema_version) || !nzchar(manifest$schema_version)) {
    add_error("schema_version must be a non-empty string")
  }
  if (is.null(manifest$artifact_root) || !nzchar(manifest$artifact_root)) {
    add_error("artifact_root must be a non-empty string")
  }
  if (is.null(manifest$runtime_classes) || length(manifest$runtime_classes) == 0L) {
    add_error("runtime_classes must be a non-empty array")
  }
  if (is.null(manifest$ci_lanes) || length(manifest$ci_lanes) == 0L) {
    add_error("ci_lanes must be a non-empty array")
  }
  if (is.null(manifest$scenarios) || length(manifest$scenarios) == 0L) {
    add_error("scenarios must be a non-empty array")
  }

  scenario_ids <- character()
  required_fields <- module_ci_required_scenario_fields()
  runtime_classes <- unlist(module_ci_coalesce(manifest$runtime_classes, character()), use.names = FALSE)
  ci_lanes <- unlist(module_ci_coalesce(manifest$ci_lanes, character()), use.names = FALSE)

  for (idx in seq_along(module_ci_coalesce(manifest$scenarios, list()))) {
    scenario <- manifest$scenarios[[idx]]
    label <- scenario$scenario_id %||% sprintf("<scenario %d>", idx)
    missing_fields <- setdiff(required_fields, names(scenario))
    if (length(missing_fields) > 0L) {
      add_error(sprintf("%s missing fields: %s", label, paste(missing_fields, collapse = ", ")))
      next
    }

    scenario_ids <- c(scenario_ids, scenario$scenario_id)
    if (!grepl("^MCI-[0-9]{3}\\.[0-9]+-[a-z0-9][a-z0-9-]*$", scenario$scenario_id)) {
      add_error(sprintf("%s has invalid scenario_id format", label))
    }
    if (!grepl("^MCI-[0-9]{3}$", scenario$ticket_id)) {
      add_error(sprintf("%s has invalid ticket_id", label))
    }
    if (!startsWith(scenario$scenario_id, paste0(scenario$item_id, "-"))) {
      add_error(sprintf("%s scenario_id must start with item_id", label))
    }
    if (!startsWith(scenario$item_id, paste0(scenario$ticket_id, "."))) {
      add_error(sprintf("%s item_id must belong to ticket_id", label))
    }
    if (!scenario$omic %in% .MODULE_CI_VALID_OMICS) {
      add_error(sprintf("%s has unsupported omic: %s", label, scenario$omic))
    }
    if (!scenario$runtime %in% runtime_classes) {
      add_error(sprintf("%s has unsupported runtime: %s", label, scenario$runtime))
    }
    if (!scenario$ci_lane %in% ci_lanes) {
      add_error(sprintf("%s has unsupported ci_lane: %s", label, scenario$ci_lane))
    }
    if (is.null(scenario$fixture$path) || !nzchar(scenario$fixture$path)) {
      add_error(sprintf("%s fixture.path must be non-empty", label))
    } else {
      fixture_path <- file.path(fixture_root, scenario$fixture$path)
      if (!file.exists(fixture_path)) {
        add_error(sprintf("%s fixture does not exist: %s", label, fixture_path))
      }
    }
    if (length(scenario$expected_artifacts) == 0L) {
      add_error(sprintf("%s expected_artifacts must be non-empty", label))
    }
    if (length(scenario$required_tests) == 0L) {
      add_error(sprintf("%s required_tests must be non-empty", label))
    }
  }

  duplicated_ids <- unique(scenario_ids[duplicated(scenario_ids)])
  if (length(duplicated_ids) > 0L) {
    add_error(sprintf("duplicate scenario_id values: %s", paste(duplicated_ids, collapse = ", ")))
  }

  if (length(errors) > 0L) {
    stop(
      sprintf(
        "Invalid module CI manifest %s:\n- %s",
        manifest_path,
        paste(errors, collapse = "\n- ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

module_ci_scenarios <- function(
    manifest = read_module_ci_manifest(),
    omic = NULL,
    module_family = NULL,
    runtime = NULL,
    ci_lane = NULL,
    ticket_id = NULL,
    scenario_id = NULL
) {
  scenarios <- manifest$scenarios
  keep <- rep(TRUE, length(scenarios))

  field_filter <- function(field, values) {
    if (is.null(values)) {
      return(rep(TRUE, length(scenarios)))
    }
    vapply(scenarios, function(scenario) scenario[[field]] %in% values, logical(1))
  }

  keep <- keep & field_filter("omic", omic)
  keep <- keep & field_filter("module_family", module_family)
  keep <- keep & field_filter("runtime", runtime)
  keep <- keep & field_filter("ci_lane", ci_lane)
  keep <- keep & field_filter("ticket_id", ticket_id)
  keep <- keep & field_filter("scenario_id", scenario_id)

  scenarios[keep]
}

module_ci_scenario_ids <- function(scenarios) {
  vapply(scenarios, `[[`, character(1), "scenario_id", USE.NAMES = FALSE)
}
