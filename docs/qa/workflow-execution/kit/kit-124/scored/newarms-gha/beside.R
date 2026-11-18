# Lays the new arms of the hosted-runner run beside the laptop arm each differs from by one choice (22 vs 8; 23 and 24
# vs 5), from aggregate.R's k124-cells.tsv tables. Reads tables only; decides nothing.
# Usage: Rscript --vanilla beside.R <newarms scored dir> <committed scored dir (holds s1-core .. s4-core)> <out.tsv>
args <- commandArgs(trailingOnly = TRUE)
rd <- function(p) utils::read.delim(p, quote = "", stringsAsFactors = FALSE, colClasses = c(arm = "character"))
new <- rd(file.path(args[[1L]], "k124-cells.tsv"))
old <- do.call(rbind, lapply(1:4, \(i) rd(file.path(args[[2L]], sprintf("s%d-core", i), "k124-cells.tsv"))))
new$source <- "runner"; old$source <- "laptop"
keep <- c("simulator", "cell", "arm", "subset", "source", "runs", "ok", "refused_designed", "crashed", "mean_k", "null05", "null05_se",
    "null01", "null01_se", "mean_fdp", "mean_fdp_se", "mean_fdp_upper", "fdp10_share", "fdp10_lower", "fdp10_upper", "rule3_limit",
    "power", "power_se", "rule1", "rule2", "rule3", "rule4", "rule5", "rule6", "rule7", "verdict", "void")
pairs <- list(c("22", "8"), c("23", "5"), c("24", "5"))
out <- list()
for (sim in paste0("S", 1:4)) for (cell in paste0("C", 0:7)) for (p in pairs) {
    n <- new[new$simulator == sim & new$cell == cell & new$arm == p[[1L]] & new$subset == "full", keep]
    o <- old[old$simulator == sim & old$cell == cell & old$arm == p[[2L]] & old$subset == "full", keep]
    if (!nrow(n)) next
    stopifnot(nrow(n) == 1L, nrow(o) == 1L)
    n$pair <- o$pair <- paste(p, collapse = " vs "); n$role <- "new"; o$role <- "reference"
    out[[length(out) + 1L]] <- rbind(n, o)
}
out <- do.call(rbind, out)
write.table(out, args[[3L]], sep = "\t", quote = FALSE, row.names = FALSE)
# arms 1 and 2: the scored figures of the two sources, largest absolute difference over the rule-bearing measures
chk <- merge(new[new$arm %in% c("1", "2"), ], old[old$arm %in% c("1", "2"), ], by = c("simulator", "cell", "arm", "subset"))
ms <- c("null05", "null01", "mean_fdp", "fdp10_share", "power", "mean_k", "anycall_share")
cat("arms 1 and 2, scored rows on both sides:", nrow(chk), "; largest absolute difference per measure:\n")
print(vapply(ms, \(m) max(abs(chk[[paste0(m, ".x")]] - chk[[paste0(m, ".y")]]), na.rm = TRUE), numeric(1L)))
