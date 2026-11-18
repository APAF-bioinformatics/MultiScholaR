# KIT-124, extra-hard cells: the option-set designs side by side, read from the final aggregate.R tables.
# Usage: Rscript option-arms.R <final dir with s1..s4/k124-cells.tsv and inventory.tsv> ; writes option-arms-extra-hard.tsv there.
# Reads only; no rule is recomputed (every status is aggregate.R's).
dir <- commandArgs(TRUE)[[1L]]
cells <- do.call(rbind, lapply(1:4, \(i) utils::read.delim(file.path(dir, paste0("s", i), "k124-cells.tsv"),
    colClasses = "character", quote = "")))
inv <- utils::read.delim(file.path(dir, "inventory.tsv"), colClasses = "character", quote = "")
design <- c("1" = "none", "2" = "ceiling (true factor)", "3" = "legacy (today)", "5" = "first-pass nulls",
    "7" = "recommended (all features, k = E)", "8" = "known/supplied controls", "11" = "arm 7 with k = 1",
    "17" = "known batch in the model")
x <- cells[cells$cell_role == "stress" & cells$subset == "full" & cells$arm %in% names(design), ]
num <- function(v) suppressWarnings(as.numeric(v))
out <- data.frame(simulator = x$simulator, cell = x$cell,
    platform = inv$platform[match(paste(x$simulator, x$cell), paste(inv$simulator, inv$cell))],
    arm = x$arm, design = unname(design[x$arm]), runs = x$runs, ok = x$ok,
    refused = num(x$refused_designed) + num(x$refused_method), crashed = x$crashed, scored = x$scored,
    mean_k = x$mean_k, type1_05 = x$null05, type1_05_se = x$null05_se, type1_01 = x$null01, type1_01_se = x$null01_se,
    mean_fdp = x$mean_fdp, mean_fdp_se = x$mean_fdp_se, fdp10_share = x$fdp10_share, fdp10_lower = x$fdp10_lower,
    fdp10_upper = x$fdp10_upper, rule3_limit = x$rule3_limit, power = x$power, power_se = x$power_se,
    rule1 = x$rule1, rule2 = x$rule2, rule3 = x$rule3, rule4 = x$rule4, rule5 = x$rule5, rule6 = x$rule6,
    ceiling_rules_1_2 = x$arm2_rule12, void = x$void,
    legacy_platform_caveat = x$arm == "3" & inv$platform[match(paste(x$simulator, x$cell), paste(inv$simulator, inv$cell))] == "runner",
    stringsAsFactors = FALSE)
all_cells <- paste0("S", 1:14)
out <- out[order(out$simulator, match(out$cell, all_cells), match(out$arm, c("2", "7", "11", "8", "5", "3", "1", "17"))), ]
utils::write.table(out, file.path(dir, "option-arms-extra-hard.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
