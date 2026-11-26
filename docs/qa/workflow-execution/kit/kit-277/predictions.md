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

## Rulings of 2026-11-26 (coordinator): new test lines, tests-only commit before the code commit
Departures recorded (step0/outcome): F2b (no real reopen) and P5 (worker path not drivable under testServer) accepted as
written, listed under "Not checked" with the lines read (`R/mod_omic_da.R` registerOmicDaDesignFormula;
`R/utils_workflow_execution.R:2842-2859` for the worker's project-copy read).
New lines and their predictions on the OLD code (R1c, new file only):
- Copy press itself (`kit277CopyClick`: Save Workflow Arguments then Copy through the Summary tab): F1/F2 (3 lanes), O4: after
  the Save "Copy error" is expected; on the old code Copy succeeds, so these expectations FAIL by value. After the re-run no
  "Copy error" and at least one note: passes on both sides (a guard).
- P4: the head rule `assertWorkflowSummaryHeadFinal` (the function the Summary Save :297, the Report :524 and Export Session
  State :622 call) must refuse with class `omic_norm_correlation_input_stale` and a text, per lane; on the old code the head is
  the final table, "allowed": FAILS by value (4 lanes). P4 prot_lfq also presses Copy: "Copy error" expected; old code fails.
  The Export Session State button itself is not pressed: its guard is the same function (listed under "Not checked").
- G4 (new, guard): an import with taxon 9606 and the name "Homo sapiens (human)" asks nothing and says "nothing changed":
  passes on both sides. (It found a defect in the code stage before any run: the organism rule asked on a name difference
  with an equal taxonomy id; the rule now decides by taxonomy id when the bus holds one.)

## SG1 (declared before sending; code commit 53e4c745, stage filter MULTISCHOLAR_KIT277_STAGE=1)
Files: test-kit277-resave-carried-through.R and test-omics-module-globals.R, one job. Expected: every S1 row (F1/F2 x3 lanes with
the Copy press, F5/F6 x3, F2b, F3, F4 x2, F7) and G2 pass; S2, S3 rows are skipped; the globals test passes (no `<<-` or
global reference added to R/mod_omic_da.R or R/mod_omic_design.R). Expected 40 job-minutes, worst case 6.67 core-hours; running
total 10.60 + R1c worst 6.67 + SG1 worst 6.67 = 23.94 of 200.00. Decision: a failing S1 row or the globals test = stop rule 3
(the one correction pass); S2/S3 do not start until SG1 holds.
