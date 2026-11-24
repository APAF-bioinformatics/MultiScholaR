# KIT-159 R1 phase pre script (r1.12): the shape of the recorded open file, as counts only.
q <- utils::read.csv(valInputPath("ref_qc"), check.names = FALSE)
utils::write.table(data.frame(item = c("rows", "columns"), value = c(nrow(q), ncol(q))),
                   file.path(out_dir, "ref-qc-shape.tsv"), sep = "\t", quote = FALSE,
                   row.names = FALSE)
