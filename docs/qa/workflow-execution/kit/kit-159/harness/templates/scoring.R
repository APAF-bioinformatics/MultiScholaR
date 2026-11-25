# KIT-NNN scoring code (section 5.3: pinned with the declaration). Ticket metrics that the harness
# does not provide, hooks and custom or published arms live here; each function the declaration
# names must be defined in a file the declaration's `code` lists.

#' A post-cell hook (section 7.3): positional arguments as revision 1.6; the named arguments are
#' passed only where the formals name them.
#' @param manager,results,summary_dir,tables As 7.3.
#' @param unit,arm,design_path,out_dir Named (revision 1.8).
#' @return Nothing; writes only the declaration's outputs.
postCell <- function(manager, results, summary_dir, tables, unit, arm, out_dir) {
    invisible(NULL)
}
