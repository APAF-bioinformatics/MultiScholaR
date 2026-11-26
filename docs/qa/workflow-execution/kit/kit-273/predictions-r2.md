# KIT-273 R2 predictions (head = the fix commit, base = tests-only 2d46e4a8; written before the run is sent)

R2 file set recounted at the branch point (the ticket's grep, no head): 45 files from the grep including the new file and
KIT-271's, plus test-gui-entry-proteomics.R = 46 (the ticket counted 45; one file entered the grep since 15d86bd4). Three
shards per side, check on both sides. Newest committed per-file minutes of the 45 existing files sum to 170.7 (ticket:
143.1); with the new file about 177 per side.

- test-kit273-limpa-ruv-delivery.R: 93/11 on base, 104/0 on head: fixed.
- test-gui-entry-proteomics.R: 92/101 on both sides (standing failure); message-changed in prot_dia_limpa's steps
  summary_copy, summary_report, summary_export only, by one added key `Study_report/design_matrix_with_ruv_factors.tab`
  each; the four file-name key lines (KIT-271's) identical on both sides; the design_matrix.tab digest does not move; every
  other record identical. If the lane's DA is not current at a Copy step the verdict is same (then reported).
- test-prot-norm-core.R 1547/307, test-workflow-saved-da.R 755/196, test-omic-da-module.R is not in this set, test-omic-qc-module.R
  not in this set: same.
- test-kit271-limpa-table-name.R, test-kit257-norm-own-files.R, test-kit244-limpa-pass-through.R, test-prot-norm-worker-figures.R,
  test-workflow-native-session.R, test-general-filemgmt-report-helpers-shared.R, test-module-ci-prot-summary.R: same.
- Every other file: same. R CMD check: identical on both sides (the function is @noRd; no Rd, no NAMESPACE).
Any other verdict stops the work (stop rule 3).
