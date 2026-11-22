# KIT-157 (branch run-k157-rc only): read a compare-results.R table and prove the rule detects a difference.
# Usage: Rscript --vanilla limpa-check.R summary <compare.tsv>
#        Rscript --vanilla limpa-check.R perturb <in.rds> <out.rds>   (one double leaf x (1 + 1e-5): must be read "outside")
arguments <- commandArgs(trailingOnly = TRUE)
if (identical(arguments[[1L]], "perturb")) {
    x <- readRDS(arguments[[2L]])
    done <- FALSE
    walk <- function(v) {
        if (done) return(v)
        if (is.list(v)) { for (i in seq_along(v)) v[[i]] <- walk(v[[i]]); return(v) }
        if (is.double(v) && length(v) && any(is.finite(v) & v != 0)) {
            i <- which(is.finite(v) & v != 0)[[1L]]
            v[[i]] <- v[[i]] * (1 + 1e-5)
            done <<- TRUE
        }
        v
    }
    x <- walk(x)
    stopifnot(done)
    saveRDS(x, arguments[[3L]])
    cat("K157CHECK perturbed one double by a relative 1e-5\n")
} else {
    t <- utils::read.delim(arguments[[2L]], stringsAsFactors = FALSE)
    doubles <- t[t$kind == "double", ]
    cat(sprintf("K157CMP leaves %d: agree %d, outside %d, differ %d; doubles %d leaves, %.0f values\n", nrow(t),
        sum(t$verdict == "agree"), sum(t$verdict == "outside"), sum(t$verdict == "differ"), nrow(doubles), sum(doubles$n, na.rm = TRUE)))
    if (nrow(doubles)) {
        a <- doubles[which.max(doubles$max_abs), ]; r <- doubles[which.max(doubles$max_rel), ]
        cat(sprintf("K157CMP largest absolute difference %.3g (%s); largest relative %.3g (%s)\n", a$max_abs, a$path, r$max_rel, r$path))
    }
    bad <- t[t$verdict != "agree", ]
    for (i in seq_len(min(nrow(bad), 40L))) {
        cat(sprintf("K157CMP %s %s %s n=%s rel=%s abs=%s outside=%s %s\n", bad$verdict[[i]], bad$path[[i]], bad$kind[[i]],
            bad$n[[i]], format(bad$max_rel[[i]], digits = 3), format(bad$max_abs[[i]], digits = 3), bad$outside[[i]], bad$note[[i]]))
    }
}
