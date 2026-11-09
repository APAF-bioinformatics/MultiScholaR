suppressMessages(library(openxlsx))
a <- commandArgs(TRUE); files <- list.files(if (length(a)) a[1] else file.path(Sys.getenv("TMPDIR"), "fixture", "da"), pattern = "_long_annot[.]xlsx$", full.names = TRUE)
used <- function(g) sum(g[, 2]); mx <- function(g) sum(g[, 6])
base <- used(gc())
stage <- function(label, expr) { gc(reset = TRUE); g0 <- used(gc()); force(expr); g <- gc(); cat(sprintf("%-14s start %6.1f  max %6.1f  after %6.1f (above base)\n", label, g0 - base, mx(g) - base, used(g) - base)) }
wb <- createWorkbook(); addWorksheet(wb, "DA_Results_Index")
for (i in seq_along(files)) {
  stage(paste0("read", i), x <- read.xlsx(files[i]))
  stage(paste0("write", i), { addWorksheet(wb, paste0("S", i)); writeData(wb, paste0("S", i), x) })
  rm(x)
}
stage("save", saveWorkbook(wb, file.path(Sys.getenv("TMPDIR"), "s.xlsx"), overwrite = TRUE))
stage("rm", rm(wb)); st <- readLines("/proc/self/status"); cat(grep("^(VmHWM|RssAnon)", st, value = TRUE), sep = "\n")
