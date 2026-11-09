# KIT-108 W2: D7's two replica rules on M2/M3's outputs. Usage: Rscript validity.R <O dir> <out.tsv>
a <- commandArgs(TRUE); O <- a[[1]]
rows <- function(d) utils::read.delim(file.path(O, d, "replica-rows.tsv"), stringsAsFactors = FALSE)
v <- function(r, phase, col = "rss_anon_kb") r[[col]][r$phase == phase]
GiB <- 1048576
out <- do.call(rbind, lapply(c("108p", "108q", "099l"), function(k) {
  reps <- if (k == "099l") "099l-i3" else paste0(k, "-i3-r", 1:2)
  w <- if (k != "099l") utils::read.delim(file.path(O, paste0(k, "-i2c-windows.tsv"))) else NULL
  gui <- if (!is.null(w)) w$r_anon_exit[w$window == "render"] - w$r_anon_entry[w$window == "copy"] else NA
  do.call(rbind, lapply(reps, function(d) {
    r <- rows(d)
    copy <- (v(r, "after_copy_raw") - v(r, "before_copy")) / GiB
    render <- (v(r, "after_render_raw") - v(r, "before_render")) / GiB
    both <- (v(r, "after_render_raw") - v(r, "before_copy")) / GiB
    data.frame(member = k, replica = d, copy_raw_rise = copy, render_raw_rise = render, replica_rise = both,
      gui_copy_entry_to_render_exit = gui, ratio_to_gui = both / gui,
      gc_used_after_copy_mb = v(r, "after_copy", "gc_used_ncells_mb") + v(r, "after_copy", "gc_used_vcells_mb"),
      gc_max_during_copy_mb = v(r, "after_copy", "gc_max_ncells_mb") + v(r, "after_copy", "gc_max_vcells_mb"),
      gc_used_before_copy_mb = v(r, "before_copy", "gc_used_ncells_mb") + v(r, "before_copy", "gc_used_vcells_mb"),
      rss_after_gc = v(r, "after_copy") / GiB, free_arena_before_render = v(r, "before_render", "free_arena_bytes") / 1024^3,
      rss_final_gc = v(r, "final_gc") / GiB, rss_after_trim = v(r, "after_malloc_trim") / GiB,
      trim_returned = (v(r, "final_gc") - v(r, "after_malloc_trim")) / GiB,
      left_above_before_copy = (v(r, "after_malloc_trim") - v(r, "before_copy")) / GiB, copy_s = v(r, "after_copy", "seconds"))
  }))
}))
num <- vapply(out, is.numeric, logical(1)); out[num] <- lapply(out[num], function(x) round(x, 4))
utils::write.table(out, a[[2]], sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 250); print(out, row.names = FALSE)
for (k in c("108p", "108q")) { x <- out[out$member == k, ]
  cat(k, "r1/r2 copy rise differ", round(100 * abs(diff(x$copy_raw_rise)) / mean(x$copy_raw_rise), 3), "%; render rise differ",
    round(100 * abs(diff(x$render_raw_rise)) / mean(x$render_raw_rise), 2), "%; replica vs GUI", round(100 * (x$ratio_to_gui - 1), 1), "%\n") }
