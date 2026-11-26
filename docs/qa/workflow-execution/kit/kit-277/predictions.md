# KIT-277 predictions (committed before R1)

R1 = the tests-only commit on the unchanged code (R/ is staging 38443af8). Binding: the commit this file is committed with.

## New file `test-kit277-resave-carried-through.R` (fails on the old code, each for its stated reason)
- F1/F2 (prot_lfq, lipid_canonical, metab_lc): `out$asked` FALSE (no question: "Design saved; nothing changed."), the DA
  tab's box keeps `~ 0 + group` (proteomics: refilled with the old formula), the DA reference stays, Copy "allowed".
- F5, F6: fail through F1's first expectation (`asked`); the notice is "Design saved; nothing changed.".
- F2b: `omicDaDesignFormula` does not exist: the read is an error text, the expectation fails; the run fits `~ 0 + group`.
- F3 (trx_screen): the notice is "Design saved; nothing changed." (fails on the text only); the DA formula guard passes.
- F4 (prot_lfq, lipid_canonical): no question; the box keeps the old formula.
- F7: no question (the comparisons-only path), the box keeps the old formula.
- G2: passes (KIT-247's "Comparisons saved").
- O1, O3, O4, O2b, O2 (prot_dia): no question, taxon stays 9606, every planted annotation stays, Copy "allowed", the DA
  and enrichment runs use the earlier table (Organism not "Synthetic test organism").
- P1 to P9 (prot_dia chain; P4 on prot_lfq, prot_tmt, prot_dia, prot_dia_limpa; P5/P6 on prot_lfq and prot_dia_limpa):
  no question, nothing applied (the earlier table stays, no GENEX1/GENEY1/GENEZ1 in the bus table), the recorded `fasta`
  input unchanged; P2 (same bytes again) passes on both sides on the notice text and the unchanged table; P7: the older
  table stays (no GENEX1) while the new file is recorded; P8 fails (`limited` is not an argument of the old gate);
  R1: a confirmed full import refused at the checkpoint leaves the folder's RDS copies in `source_dir`; R1b: no question,
  so the row fails at `out$asked`; G3 passes on both sides.
## Edited rows of `test-kit247-design-resave.R`
- W288a (prot_lfq, lipid_canonical) and W288b: fail on the old code (no question).
## Stage filter
`MULTISCHOLAR_KIT277_STAGE=1` (SG1): the S1 rows and G2 run; S2/S3 rows are skipped.

## Instrument repair 1 of R1 (declared before sending; mechanics only, no expected value changed)
R1 run 45777669520 showed, in the new file, ten errors with no result produced: (a) F4 (two rows): `kit277ExportedDesign` and
`kit277Import` were defined below their first use at file level; moved above. (b) P1/P2, P3/P9, P4 (lfq, tmt, dia), P5/P6
(lfq), R1b: the lane step id is `qc_protein_cleanup` / `qc_peptide_accession_cleanup`, so `grep("accession_cleanup")` found
none; now `grep("cleanup")`. (c) P8: on the old code the call raised an "unused argument" error inside `expect_error(class=)`,
reported as an error, not a failure; the call is captured with `tryCatch` and `expect_s3_class(condition,
"omics_design_resave_unconfirmed")` asserts the same class. The 47+93 failures of the file and the 16 of the KIT-247 file
(W288a x2 = 10, W288b = 6) were each at the stated value; they are not re-run. Re-send: the new file only.
