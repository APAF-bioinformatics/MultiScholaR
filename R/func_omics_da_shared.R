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
