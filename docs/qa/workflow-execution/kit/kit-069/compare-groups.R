# usage: Rscript compare-groups.R <before_root> <after_root>
# A module-CI test file writes one project per test case under a random id; the ids differ between runs,
# so files are grouped by their path with the id masked and every before file must match some after file
# of its group: by bytes, or at content level (RDS with the numeric `timestamp` fields and POSIX times
# removed and paths masked; xlsx sheet by sheet; text with timestamps and paths masked).
args <- commandArgs(TRUE); b <- args[1]; a <- args[2]
norm <- function(p) gsub("[0-9a-f]{10,}", "ID", gsub("([A-Za-z0-9_.-]*?-)[0-9a-f]{6,}(/|$)", "\\1ID\\2", p))
lst <- function(root) { f <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE); f <- f[!grepl("MultiScholaR\\.so$|/pkgload", f)]; split(f, norm(substring(f, nchar(root) + 2L))) }
B <- lst(b); A <- lst(a)
scrub <- function(x) {
  if (inherits(x, "POSIXt") || inherits(x, "Date")) return("<time>")
  if (is.environment(x)) return("<env>"); if (is.function(x)) return("<fn>")
  if (is.character(x)) { x <- gsub("/home/\\S+|/tmp/\\S+", "<path>", x); return(gsub("\\d{4}-\\d{2}-\\d{2}[ T]?\\d{0,2}:?\\d{0,2}:?\\d{0,2}", "<ts>", x)) }
  if (isS4(x)) { at <- attributes(x); for (n in names(at)) if (n != "class") attr(x, n) <- scrub(at[[n]]); return(x) }
  if (is.list(x)) { if (!is.null(names(x))) x <- x[names(x) != "timestamp"]; x[] <- lapply(x, scrub); return(x) }
  x
}
mtext <- function(t) gsub("\\d{4}-\\d{2}-\\d{2}[ T]?\\d{0,2}:?\\d{0,2}:?\\d{0,2}(\\.\\d+)?( [A-Z]{2,4})?", "<ts>", gsub("/home/\\S+|/tmp/\\S+", "<path>", t))
same <- function(pb, pa) {
  if (identical(unname(tools::md5sum(pb)), unname(tools::md5sum(pa)))) return("bytes")
  ext <- tolower(tools::file_ext(pb))
  if (ext == "rds") return(if (isTRUE(all.equal(scrub(readRDS(pb)), scrub(readRDS(pa))))) "content" else NA_character_)
  if (ext == "xlsx") { sb <- readxl::excel_sheets(pb); sa <- readxl::excel_sheets(pa); if (!identical(sb, sa)) return(NA_character_)
    for (s in sb) if (!isTRUE(all.equal(as.data.frame(suppressMessages(readxl::read_excel(pb, s))), as.data.frame(suppressMessages(readxl::read_excel(pa, s)))))) return(NA_character_); return("content") }
  if (ext %in% c("txt","tab","tsv","csv","json","ini","md","log")) return(if (identical(mtext(readLines(pb, warn = FALSE)), mtext(readLines(pa, warn = FALSE)))) "content" else NA_character_)
  NA_character_
}
tally <- c(bytes = 0L, content = 0L, unmatched = 0L); unmatched <- character()
for (k in names(B)) { for (pb in B[[k]]) { hits <- vapply(A[[k]] %||% character(), function(pa) same(pb, pa), character(1)); how <- if ("bytes" %in% hits) "bytes" else if ("content" %in% hits) "content" else "unmatched"; tally[how] <- tally[how] + 1L; if (how == "unmatched") unmatched <- c(unmatched, substring(pb, nchar(b) + 2L)) } }
cat(sprintf("before files %d: byte-matched %d, content-matched (timestamps/paths/ids masked) %d, unmatched %d; after-only groups %d\n", sum(lengths(B)), tally["bytes"], tally["content"], tally["unmatched"], length(setdiff(names(A), names(B)))))
for (u in unmatched) cat("  UNMATCHED", u, "\n")
