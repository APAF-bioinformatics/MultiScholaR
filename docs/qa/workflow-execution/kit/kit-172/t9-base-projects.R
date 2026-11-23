# KIT-172 T9 (ticket section 6.4, A18): the base-saved projects, written once on a runner by the base commit.
#
# Usage (run-tests scripts lane, from the package root of the commit under test):
#   Rscript docs/qa/workflow-execution/kit/kit-172/t9-base-projects.R <base sha> <out dir>
# A git worktree of the base commit is added under RUNNER_TEMP (or a temporary directory); a child Rscript, whose working
# directory is that worktree, loads it with pkgload::load_all (its own tests/testthat helpers in an environment under
# the global one, plus this commit's helper-kit172.R for fixture F172) and saves one small lipid and one small metab
# project at the base's default offset (1): F172's first 40 features and its 5 internal standards, intensity filter,
# QC finalise, ITSD, log2, cyclic loess, the analysis record. The two project roots are packed in one tar.gz and
# written base64-encoded to <out dir>/base-projects.tar.gz.b64.txt (an archive is withheld from the upload); the
# base SHA and the child's own `git rev-parse HEAD` go to <out dir>/base-projects.sha and <out dir>/status.tsv.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, grepl("^[0-9a-f]{40}$", args[[1L]]))
base <- args[[1L]]
head_root <- normalizePath(".")
out <- normalizePath(args[[2L]], mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
status <- function(step, outcome) cat(sprintf("%s\t%s\n", step, gsub("[\t\n]", " ", outcome)), file = file.path(out, "status.tsv"), append = TRUE)
cat("step\toutcome\n", file = file.path(out, "status.tsv"))

scratch <- Sys.getenv("RUNNER_TEMP", tempdir())
worktree <- file.path(scratch, "kit172-base")
if (system2("git", c("worktree", "add", "--detach", shQuote(worktree), base)) != 0L) {
  status("worktree", "mechanism: git worktree add failed")
  quit(save = "no", status = 0L)
}
projects <- file.path(scratch, "kit172-projects")
dir.create(projects, showWarnings = FALSE)
child <- file.path(scratch, "kit172-t9-child.R")
writeLines(c(
  'args <- commandArgs(trailingOnly = TRUE); head_root <- args[[1L]]; projects <- args[[2L]]',
  'suppressMessages(pkgload::load_all(".", export_all = TRUE, helpers = FALSE, attach_testthat = FALSE, quiet = TRUE))',
  'env <- new.env(parent = globalenv())',
  'suppressMessages(testthat::source_test_helpers("tests/testthat", env = env))',
  'sys.source(file.path(head_root, "tests", "testthat", "helper-kit172.R"), envir = env)',
  'cat(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), file = file.path(projects, "child-sha.txt"))',
  'local(envir = env, {',
  '  values <- kit172Values(172001L)$values[c(1:40, 516:520), ]',
  '  for (class in c("lipid", "metab")) {',
  '    root <- file.path(projects, class)',
  '    offset <- kit172DefaultOffset(class)',
  '    stopifnot(identical(offset, 1))',
  '    journey <- kit172Journey(root, kit172Object(values, class), class, da = FALSE, offset = offset)',
  '    writeWorkflowAnalysisRecord(journey$owner, file.path(root, "analysis_record.json"))',
  '    closeDIAStateAdapter(journey$owner)',
  '  }',
  '})'), child)
code <- system2("bash", c("-c", shQuote(sprintf("cd %s && Rscript %s %s %s", shQuote(worktree), shQuote(child), shQuote(head_root),
  shQuote(projects)))))
child_sha <- if (file.exists(file.path(projects, "child-sha.txt"))) readLines(file.path(projects, "child-sha.txt"), warn = FALSE) else "<none>"
status("child_sha", child_sha)
if (code != 0L || !identical(child_sha, base) || !all(file.exists(file.path(projects, c("lipid", "metab"), "analysis_record.json")))) {
  status("save", sprintf("mechanism: child exit %d, child sha %s", code, child_sha))
  quit(save = "no", status = 0L)
}
archive <- file.path(scratch, "base-projects.tar.gz")
old <- setwd(projects)
utils::tar(archive, files = c("lipid", "metab"), compression = "gzip", tar = "internal")
setwd(old)
writeLines(jsonlite::base64_enc(readBin(archive, "raw", file.size(archive))), file.path(out, "base-projects.tar.gz.b64.txt"))
writeLines(c(base, digest::digest(file = archive, algo = "sha256")), file.path(out, "base-projects.sha"))
status("save", sprintf("pass: %d bytes", file.size(archive)))
