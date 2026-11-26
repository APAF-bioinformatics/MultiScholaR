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

## R2 (declared before sending): head = the code commit, base = the tests-only tree f16edb13 (R/ = staging 38443af8)
134 files, 4 shards per side in the committed order (order4.json), file_timeouts 120 (test-kit247-design-resave.R) and 90
(the new file), check true, MULTISCHOLAR_TRX_MANIFEST_NOTE=1, blas reference.
Per file (compare.tsv verdicts):
- test-kit277-resave-carried-through.R: base fails (R1 values), head passes every row: `fixed` (base failures, head none).
  Where a head row still fails, that is a miss (stop rule 3), not re-run.
- test-kit247-design-resave.R: base fails W288a x2 and W288b (16 expectations), head passes them: `fixed`; every other
  row passes on both sides (E1, E2, E7, E12, E17, E18, E19 unchanged).
- test-kit247-import-config.R, test-kit247-resave-no-session.R: `same`. test-omics-module-globals.R: `same`.
- Every other file (including the module-CI DA and enrichment files, test-workflow-saved-da.R, test-kit212-da-current.R,
  test-kit205-enrichment-background.R, test-prot-annotation-uniprot-empty.R, test-workflow-prot-declarations.R,
  test-prot-04-design.R, test-prot-04a-design-module-characterization.R, test-prot-design-core.R,
  test-kit216-prot-session-after-enrichment.R and the 8 files that joined the grep since 7bbc344e): `same`.
- R CMD check: identical on both sides (new functions are @noRd; no Rd, no NAMESPACE change).
- No stored value moves. Any other verdict stops the work.
Budget line: spent 10.60 + R1c (about 3) + SG1 (about 3) = about 17 before R2; R2 expected 58.04, worst case 92.00; total
worst case under 200.00.

## Diagnostic for O4's miss (declared before sending; outside the correction pass; nothing is judged by it)
Purpose: find why, on the head, O4's Copy gate says "refused" and the Copy press raises the "files unmade" text after the organism
Save, the DA run and the enrichment run (R2 45804909404). Tests-only commit: `kit277Diag` (log lines `KIT277DIAG <point> <name> = <value>`
written to a file and printed by the O4 test when MULTISCHOLAR_KIT277_DIAG is set) and a stage value "2" (S2 rows only). No expectation
changes. Points: A after the organism Save; B after the first Save Workflow Arguments + Copy; C after the refused enrichment press;
D after the DA run; E after the enrichment run; F after the second Save + Copy. Logged: manager and snapshot current state, history with
each state's entity id, context head entity, store identity, da_results entity and its manifest's source entity/owner/parameters,
enrichment marker, workflowDaResultStatus and the three tests (owner, membership, workflowDaTablesDiffer with the component ids),
the head's args$daAnalysisUI, candidate against made against unmade DA files, the publication marker, every notification, the Copy
gate's class and full message. One send: base (tests tree with the diagnostic on the staging R/) and head (the same tests with the
code), env MULTISCHOLAR_KIT277_STAGE=2, MULTISCHOLAR_KIT277_DIAG=1; file timeout 25 minutes; worst case 4.67 core-hours of the 5.0
authorised; running total about 71 + 4.67 of 200.00. Hand back the log's reading; no repair.

## Diagnostic 2 for O4's miss (declared before sending; outside the correction pass; coordinator-allowed; nothing judged by it)
Purpose: the full text and class of every notification and of any condition raised in the DA press at point D, and the DA handler's
console output for that press. Added (tests only, only when MULTISCHOLAR_KIT277_DIAG is set): `notes_full` at every point (no `str`);
at D a tracing wrapper of `runDIAReferenceDaAnalysis` (the in-process DA call at R/mod_prot_da_handlers_run_analysis.R:~138) that logs
the formula the request carries, the annotation argument (rows and columns, or class) and the taxon, and on error the condition's class,
call and message before re-raising; the press's console output (`cat()` lines) captured into the diagnostic file; the bus annotation's
rows, id column, organisms and the removed flag after the press; at E the enrichment press's console output. No expectation changes. One
send, base (diagnostic file on staging R) and head (same file, code), env STAGE=2 and DIAG=1, file timeout 25 minutes; worst case 4.0
core-hours (3.9 expected ceiling: 2 jobs x 35 min = 70 min = 4.67 would exceed, so the file timeout is lowered to 20 min: 2 x 30 = 60 min
= 4.0); running total about 74 + 4.0 of 200.00. Hand back: the error text and class, the code line, the smallest repair as text.

## Diagnostic 3 for O4's miss (declared before sending; outside the correction pass; coordinator-allowed; nothing judged by it)
Purpose: read why `export_pending` is set when the DA press runs after an organism Save (diagnostic 2: `assertDIAExportReady` refuses
with saved_entity_outputs_pending). Tests-only (only when MULTISCHOLAR_KIT277_DIAG is set): `kit277ExportState` (export_pending null or
not with status/generation/request id; head entity; `owned_export$generation` and whether it equals the head; the runner the lane uses
(`workflowSessionApplication` null or its class); the archive's member list, with the DA members); `kit277CopyClick` with a `diag`
argument that, around the Save Workflow Arguments press, logs the state before and after, the Save press's OWN notifications in full
(notes cleared before it, read before the Copy press), a trace of `refreshDIANativeSummaryArchive` (called, early return expected or
not, error class/call/message, export_pending after) and of `runDIANativeSessionExport` (called, publication status, error class/call/
message); points B3 (first Save after the organism Save), D3_pre_da_press, F3 (second Save); and a diagnostic-only test (no expectation)
for the formula Save lane prot_lfq: export state before and after its Save and its Save Workflow Arguments, one block for contrast. No
expectation is changed. One send, base and head, STAGE=2 and DIAG=1, file timeout 20 minutes; worst case 3.0 core-hours (2 x 30 min =
60 min = 4.0 exceeds, so the timeout is 15 minutes: 2 x 25 = 50 min = 3.33... lowered to a 12-minute file timeout: 2 x 22 = 44 min =
2.93); running total about 77 + 2.93 of 200.00.

## THE ONE CORRECTION PASS (coordinator's ruling 2026-11-27; AUDIT-KIT-277-plan section "Miss 1: cause and repair")
Tests-only commit 1: (i) ruling K: P4 prot_lfq's `expect_false(out$record_da)` is replaced by the user-level guarantee: the digest of
`source_dir/analysis_record.json` is unchanged across the Save Workflow Arguments press after the FASTA cut (the press is refused by the head
rule and writes no record); (ii) new row O5 (S2): the session exported at the DA head (the lane's session_export step), organism Save, Save
Workflow Arguments: `export_pending` NULL; the Save's notices hold no error; the owned generation equals the head; the refreshed metadata
holds no `uniprot_dat_cln`; the members of `collectDIANativeGuiFiles` include no `source/uniprot_dat_cln.RDS` and no `da/`; the
annotation file is still absent. `kit277CopyClick` also returns the Save press's notices as an attribute and clears notes before the Save.
(iii) The diagnostic blocks (kit277Diag, kit277ExportState, the traces, the diagnostic-only test) STAY in the file, gated by
MULTISCHOLAR_KIT277_DIAG, which is unset in every counted run (R1, R2): declared.
Code commit 2: in `refreshDIANativeSummaryArchive` (R/utils_workflow_dia_sessions.R), before the export call, take `aa_seq_tbl_final` and
`uniprot_dat_cln` from the CURRENT bus into the metadata, and in its `prepare_outputs` write those two tables as the user's export does
(R/mod_prot_norm_session_helpers.R:256-261); nothing else. Option A is rejected. The DA readiness check, the pending flag's lifecycle and
`setOmicsDaBack` are not touched.
Predictions, R1 (the tests-only tree; new file and test-kit247-design-resave.R, old code): the new file: O5's lines fail by value on
the old code (the earlier archive still declares the annotation and the DA folder: members include `source/uniprot_dat_cln.RDS` and
`da/`, annotation_in_metadata TRUE, source file present); P4 prot_lfq's replaced line fails by value (the old code applies no cut, so
the Save writes a new analysis_record.json); O4's two lines PASS on the old code as in R2; every other row as in R1/R1c (fail by value,
except G2, G3, G4, P2 parts, P4 head-rule for tmt/dia/limpa which pass on both sides). The KIT-247 file: 170 passed / 16 failed as in R1.
Predictions, R2 (head = the code commit, base = the tests-only tree; the 134 files plus test-gui-entry-proteomics.R, 4 shards per side
from order4.json with the gui-entry file in the lightest shard, check both, file timeouts 120 (KIT-247 file, gui-entry file) and 90 (new
file), MULTISCHOLAR_TRX_MANIFEST_NOTE=1): the new file passes every row on the head (O4 and O5 included); the KIT-247 file `fixed`
170/16 -> 186/0; test-gui-entry-proteomics.R `same` (92 passed / 101 as on staging; the refresh runs in every proteomics lane's Save
path: a stored GUI-entry leaf recording source/aa_seq_tbl_final.RDS or source/uniprot_dat_cln.RDS moving is a RESULT and stops the work);
test-kit216-prot-session-after-enrichment.R, test-kit206-prot-norm-export-guard.R, test-workflow-native-session.R and every other file
`same`; check identical. A miss anywhere in R2: stop, the correction pass is spent.
Budget line: spent about 79.5 (R1..R2, three diagnostics, incl. diagnostic 3); R1 expected 7.3 (worst 10), R2 expected 56 (worst 70)
-> at most 79.5 + 10 + 70 = 159.5 of 200.00.

## R2 of the correction pass: running-ceiling line (2026-11-27, before the request is pushed)
Used so far from runs.tsv (core-hours = job-minutes x 4 / 60): 7.73 + 2.87 + 3.0 + 1.9 + 55 + 2.5 + 2.5 + 2.45 + 8.6 = 86.55 of 200.00
(R1 of the correction pass 45837030184: about 129 job-minutes = 8.6; it was read as predicted). R2: expected about 56 (the earlier R2 was
about 806 job-minutes = 53.7, plus test-gui-entry-proteomics.R 14.1 minutes a side = 1.9, plus the 120-minute timeouts only if a file stalls),
worst case 92.00 (as declared for the first R2: 8 test jobs and 2 check jobs, 4 cores each). 86.55 + 92.00 = 178.55 <= 200.00.
Request: head 174a62db9 (R/ = the correction code; tests identical to the base), base dfffaa0b (tests-only tree, R/ = staging), 135 files,
4 shards per side, round-robin order built from order4.json: test-gui-entry-proteomics.R joins the lightest shard of the earlier R2 (shard 1,
69.8 minutes); to keep round-robin sizes (34, 34, 34, 33) the lightest file of that shard (test-kit247-resave-no-session.R, 4.7 minutes) moves to
the former shard 4; shard loads after the move about 84, 95, 93, 99 minutes. File assignment does not change any result (one process per file).
file_timeouts: 120 for test-kit247-design-resave.R and test-gui-entry-proteomics.R, 90 for the new file; default 45; env
MULTISCHOLAR_TRX_MANIFEST_NOTE=1 only (no MULTISCHOLAR_KIT277_DIAG, no stage variable); check true.
