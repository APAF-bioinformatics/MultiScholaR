.MODULE_CI_FIXTURE_CLASSES <- c(
  "happy_path",
  "missingness",
  "duplicates",
  "invalid_design",
  "multi_factor",
  "bad_names",
  "nonfinite",
  "small_n",
  "large_enough_for_plots",
  "multi_assay_mixed"
)

.MODULE_CI_FIXTURE_CATEGORIES <- c("normal", "malformed", "boundary", "multi_assay", "branch_specific")
.MODULE_CI_FIXTURE_LANES <- c("push", "nightly")

module_ci_expected_import_formats <- function() {
  list(
    proteomics = c("dia", "dda", "tmt", "dia_limpa"),
    metabolomics = c("lc", "gc"),
    lipidomics = c("lc", "gc")
  )
}

module_ci_fixture_pack_path <- function() {
  file.path(module_ci_fixture_root(), "fixture_packs.json")
}

module_ci_as_character <- function(x) {
  as.character(unlist(module_ci_coalesce(x, character()), use.names = FALSE))
}

module_ci_pack_path <- function(pack, field, fixture_root = module_ci_fixture_root()) {
  value <- pack[[field]]
  if (is.null(value) || !nzchar(value)) {
    stop(sprintf("Pack %s is missing %s", pack$pack_id %||% "<unknown>", field), call. = FALSE)
  }
  file.path(fixture_root, value)
}

read_module_ci_fixture_packs <- function(
    pack_path = module_ci_fixture_pack_path(),
    fixture_root = dirname(pack_path),
    validate = TRUE
) {
  if (!file.exists(pack_path)) {
    stop(sprintf("module CI fixture pack catalog not found: %s", pack_path), call. = FALSE)
  }
  catalog <- jsonlite::read_json(pack_path, simplifyVector = FALSE)
  if (isTRUE(validate)) {
    validate_module_ci_fixture_packs(catalog, fixture_root = fixture_root, pack_path = pack_path)
  }
  packs <- catalog$packs
  names(packs) <- vapply(packs, `[[`, character(1), "pack_id")
  catalog$packs <- packs
  catalog
}

validate_module_ci_fixture_packs <- function(
    catalog,
    fixture_root = module_ci_fixture_root(),
    pack_path = module_ci_fixture_pack_path()
) {
  errors <- character()

  add_error <- function(message) {
    errors <<- c(errors, message)
  }

  required_catalog_fields <- c(
    "schema_version",
    "ticket_id",
    "generated_by",
    "deterministic_seed",
    "fixture_classes",
    "import_formats",
    "packs"
  )
  missing_catalog <- setdiff(required_catalog_fields, names(catalog))
  if (length(missing_catalog) > 0L) {
    add_error(sprintf("catalog missing fields: %s", paste(missing_catalog, collapse = ", ")))
  }
  if (!identical(catalog$ticket_id, "MCI-002")) {
    add_error("catalog ticket_id must be MCI-002")
  }

  class_names <- vapply(
    module_ci_coalesce(catalog$fixture_classes, list()),
    function(contract) contract$fixture_class %||% "",
    character(1)
  )
  missing_classes <- setdiff(.MODULE_CI_FIXTURE_CLASSES, class_names)
  extra_classes <- setdiff(class_names, .MODULE_CI_FIXTURE_CLASSES)
  if (length(missing_classes) > 0L || length(extra_classes) > 0L) {
    add_error(sprintf(
      "fixture class contract mismatch; missing=[%s], extra=[%s]",
      paste(missing_classes, collapse = ", "),
      paste(extra_classes, collapse = ", ")
    ))
  }

  expected_formats <- module_ci_expected_import_formats()
  for (omic in names(expected_formats)) {
    actual_formats <- module_ci_as_character(catalog$import_formats[[omic]])
    if (!setequal(actual_formats, expected_formats[[omic]])) {
      add_error(sprintf(
        "import format mismatch for %s; actual=[%s], expected=[%s]",
        omic,
        paste(actual_formats, collapse = ", "),
        paste(expected_formats[[omic]], collapse = ", ")
      ))
    }
  }

  packs <- module_ci_coalesce(catalog$packs, list())
  if (length(packs) == 0L) {
    add_error("packs must be a non-empty array")
  }
  pack_ids <- character()
  required_pack_fields <- c(
    "pack_id",
    "ticket_id",
    "item_ids",
    "omic",
    "import_format",
    "fixture_class",
    "case_category",
    "ci_lane",
    "push_safe",
    "data_path",
    "design_path",
    "oracle_path",
    "generator",
    "estimated_rows",
    "estimated_columns"
  )

  for (idx in seq_along(packs)) {
    pack <- packs[[idx]]
    label <- pack$pack_id %||% sprintf("<pack %d>", idx)
    missing_pack <- setdiff(required_pack_fields, names(pack))
    if (length(missing_pack) > 0L) {
      add_error(sprintf("%s missing fields: %s", label, paste(missing_pack, collapse = ", ")))
      next
    }

    pack_ids <- c(pack_ids, pack$pack_id)
    if (!grepl("^MCI-002-[a-z]+-[a-z0-9-]+-[a-z0-9-]+$", pack$pack_id)) {
      add_error(sprintf("%s has invalid pack_id format", label))
    }
    if (!identical(pack$ticket_id, "MCI-002")) {
      add_error(sprintf("%s ticket_id must be MCI-002", label))
    }
    if (!pack$omic %in% names(expected_formats)) {
      add_error(sprintf("%s has unsupported omic: %s", label, pack$omic))
    } else if (!pack$import_format %in% expected_formats[[pack$omic]]) {
      add_error(sprintf("%s has unsupported import_format: %s", label, pack$import_format))
    }
    if (!pack$fixture_class %in% .MODULE_CI_FIXTURE_CLASSES) {
      add_error(sprintf("%s has unsupported fixture_class: %s", label, pack$fixture_class))
    }
    if (!pack$case_category %in% .MODULE_CI_FIXTURE_CATEGORIES) {
      add_error(sprintf("%s has unsupported case_category: %s", label, pack$case_category))
    }
    if (!pack$ci_lane %in% .MODULE_CI_FIXTURE_LANES) {
      add_error(sprintf("%s has unsupported ci_lane: %s", label, pack$ci_lane))
    }
    if (!is.logical(pack$push_safe) || length(pack$push_safe) != 1L) {
      add_error(sprintf("%s push_safe must be a scalar logical", label))
    }
    if (identical(pack$push_safe, FALSE) && !identical(pack$ci_lane, "nightly")) {
      add_error(sprintf("%s non-push-safe fixture must be in nightly lane", label))
    }
    if (identical(pack$fixture_class, "large_enough_for_plots") && isTRUE(pack$push_safe)) {
      add_error(sprintf("%s large plot fixtures must not be push-safe", label))
    }
    for (field in c("data_path", "design_path", "oracle_path")) {
      resolved_path <- file.path(fixture_root, pack[[field]])
      if (!file.exists(resolved_path)) {
        add_error(sprintf("%s %s does not exist: %s", label, field, resolved_path))
      }
    }
  }

  duplicated_pack_ids <- unique(pack_ids[duplicated(pack_ids)])
  if (length(duplicated_pack_ids) > 0L) {
    add_error(sprintf("duplicate pack_id values: %s", paste(duplicated_pack_ids, collapse = ", ")))
  }

  for (omic in names(expected_formats)) {
    for (import_format in expected_formats[[omic]]) {
      matching <- Filter(
        function(pack) identical(pack$omic, omic) && identical(pack$import_format, import_format),
        packs
      )
      actual_classes <- vapply(matching, `[[`, character(1), "fixture_class")
      if (!setequal(actual_classes, .MODULE_CI_FIXTURE_CLASSES)) {
        add_error(sprintf(
          "fixture class coverage mismatch for %s/%s; actual=[%s]",
          omic,
          import_format,
          paste(sort(actual_classes), collapse = ", ")
        ))
      }
    }
  }

  if (length(errors) > 0L) {
    stop(
      sprintf(
        "Invalid module CI fixture pack catalog %s:\n- %s",
        pack_path,
        paste(errors, collapse = "\n- ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

module_ci_fixture_packs <- function(
    catalog = read_module_ci_fixture_packs(),
    omic = NULL,
    import_format = NULL,
    fixture_class = NULL,
    case_category = NULL,
    ci_lane = NULL,
    push_safe = NULL
) {
  packs <- catalog$packs
  keep <- rep(TRUE, length(packs))

  field_filter <- function(field, values) {
    if (is.null(values)) {
      return(rep(TRUE, length(packs)))
    }
    vapply(packs, function(pack) pack[[field]] %in% values, logical(1))
  }

  keep <- keep & field_filter("omic", omic)
  keep <- keep & field_filter("import_format", import_format)
  keep <- keep & field_filter("fixture_class", fixture_class)
  keep <- keep & field_filter("case_category", case_category)
  keep <- keep & field_filter("ci_lane", ci_lane)
  if (!is.null(push_safe)) {
    keep <- keep & vapply(packs, function(pack) identical(pack$push_safe, push_safe), logical(1))
  }

  packs[keep]
}

module_ci_pack_ids <- function(packs) {
  vapply(packs, `[[`, character(1), "pack_id", USE.NAMES = FALSE)
}

module_ci_read_pack_table <- function(pack, field = "data_path", fixture_root = module_ci_fixture_root()) {
  utils::read.delim(
    module_ci_pack_path(pack, field, fixture_root = fixture_root),
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

module_ci_read_pack_oracle <- function(pack, fixture_root = module_ci_fixture_root()) {
  jsonlite::read_json(
    module_ci_pack_path(pack, "oracle_path", fixture_root = fixture_root),
    simplifyVector = FALSE
  )
}
