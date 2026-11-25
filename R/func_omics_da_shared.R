# ----------------------------------------------------------------------------
# getCountsTable
# ----------------------------------------------------------------------------
#' Return the quantitative table owned by a supported omics object
#'
#' @param obj A lipidomics, metabolomics, or protein quantitative-data object.
#'
#' @return The object's assay list or protein quantitative table.
#' @keywords internal
#' @noRd
getCountsTable <- function(obj) {
  if (inherits(obj, "MetaboliteAssayData")) {
    return(.getCountsTableMetabolomics(obj))
  }

  .getCountsTableLipidomics(obj)
}

# ----------------------------------------------------------------------------
# The Benjamini-Hochberg fallback of the q-value column (KIT-268, owner's decisions VAL-312 and VAL-314)
# ----------------------------------------------------------------------------
#' The sentence that says what the fallback values are, for the analysis record (the six report templates carry the
#' same text: their `ReportFdrSentence` is a copy per template). Basis: in qvalue 2.44.0 q = pi0 x min(1, u), so the
#' Benjamini-Hochberg value is q / pi0 at every p-value; for these comparisons pi0 could not be estimated.
#' @noRd
daQvalueFallbackSentence <- function() {
  paste("For these comparisons the q-value column holds Benjamini-Hochberg adjusted p-values, which equal Storey's q-values",
    "with the proportion of unchanged features (pi0) set to 1; because pi0 could not be estimated for these comparisons,",
    "the values are larger than Storey's q-values would be by a factor of 1/pi0 of unknown size (1 or more), so they are",
    "more conservative than, and not on the same scale as, the q-values of the other comparisons.")
}

#' The (assay, key) pairs of a tab's `qvalue_warnings` (assay -> comparison -> TRUE); a flat list (comparison -> TRUE) has
#' no assay. NULL when nothing fell back.
#' @noRd
daQvalueFallbackKeys <- function(qvalue_warnings) {
  rows <- lapply(names(qvalue_warnings), function(assay) {
    item <- qvalue_warnings[[assay]]
    if (!is.list(item)) return(data.frame(assay = NA_character_, key = assay, stringsAsFactors = FALSE))
    if (!length(item)) return(NULL)
    data.frame(assay = rep(assay, length(item)), key = names(item), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' The friendly name of each fallback key: an exact match on the contrasts table's `contrasts`, else the text before the
#' first "=" (the package's contrast form is `<friendly>=<contrast>`), else the key itself.
#' @noRd
daQvalueFallbackNames <- function(keys, contrasts) {
  keys <- as.character(keys)
  known <- is.data.frame(contrasts) && all(c("contrasts", "friendly_names") %in% names(contrasts))
  vapply(keys, function(key) {
    hit <- if (known) match(key, contrasts$contrasts) else NA_integer_
    if (!is.na(hit)) return(as.character(contrasts$friendly_names[[hit]]))
    if (grepl("=", key, fixed = TRUE)) sub("=.*$", "", key) else key
  }, character(1), USE.NAMES = FALSE)
}

#' The notice's text: the comparisons, grouped by assay ("A_vs_B, C_vs_D in LCMS_Neg; ..."), after what happened to them.
#' @noRd
daQvalueFallbackText <- function(comparisons, assays) {
  comparisons <- as.character(comparisons)
  assays <- rep_len(ifelse(is.na(assays), "", as.character(assays)), length(comparisons))
  groups <- split(comparisons, factor(assays, levels = unique(assays)))
  parts <- vapply(names(groups), function(assay) {
    named <- paste(unique(groups[[assay]]), collapse = ", ")
    if (nzchar(assay)) paste0(named, " in ", assay) else named
  }, character(1))
  paste0("Storey's q-values could not be computed for these comparisons, so Benjamini-Hochberg adjusted p-values are ",
    "shown in the q-value column instead: ", paste(parts, collapse = "; "), ".")
}

#' Show the notice on the tab: persistent, like the other warnings of the DA tabs.
#' @noRd
notifyDaQvalueFallback <- function(comparisons, assays) {
  shiny::showNotification(daQvalueFallbackText(comparisons, assays), type = "warning", duration = NULL)
}
