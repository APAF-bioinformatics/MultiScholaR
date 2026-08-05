#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(yaml))

source(file.path("tools", "refactor", "fidelity_inventory.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--manifest") {
  stop("Usage: convert_methods_to_helpers.R --manifest <path>", call. = FALSE)
}

manifest <- yaml::read_yaml(args[[2L]])
conversions <- manifest$method_helper_conversions
if (!length(conversions)) {
  stop("Manifest has no method_helper_conversions", call. = FALSE)
}

method_definition <- function(expr) {
  parts <- as.list(expr)[-1L]
  part_names <- names(parts) %||% rep("", length(parts))
  index <- match("definition", part_names)
  if (is.na(index) || !is.call(parts[[index]]) ||
      !identical(normalize_call_name(parts[[index]]), "function")) {
    stop("setMethod expression has no function-valued definition", call. = FALSE)
  }
  parts[[index]]
}

find_conversion <- function(inventory, conversion) {
  matches <- Filter(function(record) {
    call <- record$call_kind
    value <- record$call_value
    signature <- record$call_signature %||% ""
    identical(call, "setMethod") &&
      identical(value, conversion$method) &&
      identical(signature, conversion$signature)
  }, inventory)

  if (length(matches) != 1L) {
    stop(
      sprintf(
        "%s: expected one %s::%s method, found %d",
        conversion$source,
        conversion$method,
        conversion$signature,
        length(matches)
      ),
      call. = FALSE
    )
  }
  matches[[1L]]
}

rewrite_expression <- function(lines, record, helper_name) {
  start <- line_from_srcref(record$srcref, 1L)
  end <- line_from_srcref(record$srcref, 3L)
  block <- lines[start:end]

  definition_line <- grep("^[[:space:]]*definition[[:space:]]*=", block)[1L]
  if (is.na(definition_line)) {
    stop("Could not locate definition argument for ", helper_name, call. = FALSE)
  }

  function_line <- definition_line
  if (!grepl("function[[:space:]]*\\(", block[[function_line]])) {
    candidates <- grep("function[[:space:]]*\\(", block)
    candidates <- candidates[candidates > definition_line]
    if (!length(candidates)) {
      stop("Could not locate method function for ", helper_name, call. = FALSE)
    }
    function_line <- candidates[[1L]]
  }

  function_text <- block[[function_line]]
  if (function_line == definition_line) {
    function_text <- sub(
      "^[[:space:]]*definition[[:space:]]*=[[:space:]]*",
      paste0(helper_name, " <- "),
      function_text
    )
  } else {
    function_text <- sub(
      "^[[:space:]]*(function[[:space:]]*\\()",
      paste0(helper_name, " <- \\1"),
      function_text
    )
  }

  function_tail <- if (function_line < length(block) - 1L) {
    block[seq.int(function_line + 1L, length(block) - 1L)]
  } else {
    character()
  }
  replacement <- c(function_text, function_tail)

  roxygen_index <- start - 1L
  while (roxygen_index >= 1L && grepl("^#'", lines[[roxygen_index]])) {
    if (identical(trimws(lines[[roxygen_index]]), "#' @export")) {
      lines[[roxygen_index]] <- "#' @keywords internal"
      replacement <- c("#' @noRd", replacement)
    }
    roxygen_index <- roxygen_index - 1L
  }

  before <- if (start > 1L) lines[seq_len(start - 1L)] else character()
  after <- if (end < length(lines)) lines[seq.int(end + 1L, length(lines))] else character()
  c(before, replacement, after)
}

grouped <- split(conversions, vapply(conversions, `[[`, character(1L), "source"))

for (source_path in names(grouped)) {
  file_conversions <- grouped[[source_path]]
  inventory <- build_inventory(source_path)
  records <- lapply(file_conversions, function(conversion) {
    record <- find_conversion(inventory, conversion)
    list(
      conversion = conversion,
      record = record,
      function_hash = hash_text(ast_fingerprint(method_definition(record$expr)))
    )
  })

  starts <- vapply(records, function(item) line_from_srcref(item$record$srcref, 1L), integer(1L))
  records <- records[order(starts, decreasing = TRUE)]
  lines <- readLines(source_path, warn = FALSE)

  for (item in records) {
    lines <- rewrite_expression(lines, item$record, item$conversion$helper)
  }
  writeLines(lines, source_path)

  rewritten <- build_inventory(source_path)
  for (item in records) {
    helper <- item$conversion$helper
    helper_records <- Filter(function(record) identical(record$symbol, helper), rewritten)
    if (length(helper_records) != 1L) {
      stop(source_path, ": helper conversion failed for ", helper, call. = FALSE)
    }
    helper_expr <- helper_records[[1L]]$expr[[3L]]
    helper_hash <- hash_text(ast_fingerprint(helper_expr))
    if (!identical(item$function_hash, helper_hash)) {
      stop(source_path, ": function AST drift after converting ", helper, call. = FALSE)
    }
  }
}

cat(sprintf("Converted %d S4 methods to AST-identical helpers.\n", length(conversions)))
