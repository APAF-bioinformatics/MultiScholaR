# ----------------------------------------------------------------------------
# createOutputDir
# ----------------------------------------------------------------------------
##################################################################################################################
# =====================================================================================================
#' @export
createOutputDir <- function(output_dir, no_backup) {
    if (output_dir == "") {
        logerror("output_dir is an empty string")
        q()
    }
    if (dir.exists(output_dir)) {
        if (no_backup) {
            unlink(output_dir, recursive = TRUE)
        } else {
            backup_name <- paste(output_dir, "_prev", sep = "")
            if (dir.exists(backup_name)) {
                unlink(backup_name, recursive = TRUE)
            }
            system(paste("mv", output_dir, backup_name))
        }
    }
    dir.create(output_dir, recursive = TRUE)
}

