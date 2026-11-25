# KIT-224 control of the 12 `provenance` fields (ticket section 6; the mechanism of KIT-190 norm-dump.R, copied): test-omic-norm-module.R run on base and on head with its dump switch
# (MULTISCHOLAR_KIT105_W1_DUMP, the file's own, :38-42), and the two dumps compared field by field. Evidence code: no
# package code. Each side runs in a fresh Rscript (pkgload::load_all(<tree>, export_all = TRUE, helpers = TRUE)).
#
# Usage, from the package root of the head commit on a run-tests runner:
#   Rscript docs/qa/workflow-execution/kit/kit-224/norm-dump.R <base commit> <out dir>
# Written: base.R and head.R (the dumps), norm-module-control.tsv (path, base, head, same), sides.tsv (exit status).
args <- commandArgs(TRUE)
if (length(args) && identical(args[[1L]], "--side")) {
  tree <- args[[2L]]
  Sys.setenv(MULTISCHOLAR_KIT105_W1_DUMP = args[[3L]])
  suppressPackageStartupMessages(pkgload::load_all(tree, export_all = TRUE, helpers = TRUE, quiet = TRUE))
  testthat::test_file(file.path(tree, "tests", "testthat", "test-omic-norm-module.R"), reporter = "summary",
    stop_on_failure = FALSE)
  quit(status = 0L)
}
stopifnot(length(args) == 2L, grepl("^[0-9a-f]{40}$", args[[1L]]))
out <- args[[2L]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))
head_root <- normalizePath(".")
base_root <- file.path(tempdir(), "kit224-dump-base")
invisible(system2("git", c("worktree", "add", "--detach", shQuote(base_root), args[[1L]]), stdout = TRUE, stderr = TRUE))
stopifnot(dir.exists(file.path(base_root, "R")))
rscript <- file.path(R.home("bin"), "Rscript")
status <- c(base = NA_integer_, head = NA_integer_)
for (side in names(status)) {
  tree <- if (side == "base") base_root else head_root
  dump <- file.path(out, paste0(side, ".R"))
  unlink(dump)
  status[[side]] <- system2(rscript, c(shQuote(script), "--side", shQuote(tree), shQuote(dump)),
    stdout = file.path(out, paste0(side, ".log")), stderr = file.path(out, paste0(side, ".log")))
}
utils::write.table(data.frame(side = names(status), exit = status, commit = c(args[[1L]],
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE))), file.path(out, "sides.tsv"), sep = "\t", quote = FALSE,
  row.names = FALSE)
read <- function(side) {
  env <- new.env()
  sys.source(file.path(out, paste0(side, ".R")), envir = env)
  as.list(env)
}
flatten <- function(x, prefix = "") {
  if (is.list(x)) {
    nm <- names(x)
    if (is.null(nm)) nm <- as.character(seq_along(x))
    do.call(c, Map(function(v, n) flatten(v, if (nzchar(prefix)) paste(prefix, n, sep = "$") else n), x, nm))
  } else stats::setNames(paste(format(x, digits = 17), collapse = ","), prefix)
}
b <- flatten(read("base"))
h <- flatten(read("head"))
paths <- union(names(b), names(h))
rows <- data.frame(path = paths, base = unname(b[paths]), head = unname(h[paths]), stringsAsFactors = FALSE)
rows$same <- !is.na(rows$base) & !is.na(rows$head) & rows$base == rows$head
utils::write.table(rows, file.path(out, "norm-module-control.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("KIT-224 norm-module control: %d fields, %d differ: %s\n", nrow(rows), sum(!rows$same),
  paste(rows$path[!rows$same], collapse = "; ")))
