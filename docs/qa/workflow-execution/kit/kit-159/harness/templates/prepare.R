# KIT-NNN preparation (section 7.8 step 5). Sourced by val-run.R in an isolated environment after
# the inputs are fetched and verified; val-out/ does not exist yet. Read only declared roles
# (msrval::valInputPath); write only under work_dir.

#' The file sets per unit, the design and the contrasts (and optional declared result tables)
#' @param data_dir The verified inputs.
#' @param work_dir The preparation's work folder.
#' @return list(units = <named list: unit id -> named list of file roles to paths>,
#'   features = <data.frame(assay, feature_id): every feature of the declared input role, the key
#'   universe's source (4.1)>, design = <path>, contrasts = <list of group pairs>,
#'   result_tables = <optional named list>)
prepareInputs <- function(data_dir, work_dir) {
    stop("fill in")
}

#' The truth tables: key, class, expected, truth_status, reason (section 4.1)
#' @param data_dir The verified inputs.
#' @return Named list of data frames, one per truth_id.
buildTruth <- function(data_dir) {
    stop("fill in")
}

#' The scoring key of a feature of an arm, with its status and reason (section 4.1)
#' @param feature_id Feature identifiers.
#' @param arm Arm id (`NA` when the harness maps `features` into the key universe).
#' @return data.frame(key, status in ok/unmapped/ambiguous, reason).
featureKey <- function(feature_id, arm) {
    data.frame(key = feature_id, status = "ok", reason = NA_character_)
}
