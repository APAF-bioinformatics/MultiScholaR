# KIT-222 outcome (ticket r1.3; owner's OD-1 to OD-8 of 2026-11-25)

State: stopped after the audit-repair run (run 45574235919), as the coordinator ruled: one result not as predicted.
The test file passes in full on head (G1 included), every file is `same` and check is identical, but C1 finds 48
four-sample `manual` sequences, base merge-free, whose labels equal base and whose replicate column differs in type only
(the audit predicted C1 identical in all 52,384 base merge-free sequences; observed 52,336). The coordinator accepted the
48 as a stated limit (2026-11-25; A3). Ready for the merge train; the full suite (A5) is satisfied by the first
full-suite checkpoint on staging that contains the merge, and the ticket is not closed before it.

Recorded limits:
- E13 (I03): tested at the mechanism only; the smallest full test is named under "Deviations".
- I07: the line in `integrationReadInputs` sits before the existing validity check, so a malformed manifest
  (`biological_sample_column` of length other than 1) on a samples table with `tech_reps` gets a base-R error instead of
  the classed `integration_input_invalid`; moving the two lines below the validity block fixes it at the next touch.
  The refusal's positive case (a pre-fix input with D1) has no test. Development builds only.
- I09 (F5): NOT done; with the owner. After "Assign Technical Replicates" the selection stays filled, so adding the next
  sample's runs without removing the first ones marks two biological samples as one, and nothing refuses it. A fix
  needs the builder's selection observer changed as well (about five lines); its own small ticket if the owner says yes.

## Commits (branch kit-222, from origin/v0.9-staging cf8af245)

| Commit | What |
|---|---|
| e9309b9d | tests first: E1 to E16b, E3/E3b corpus, E7 (every builder and the applet), consumers, G1, G4, G5 |
| bcdd8bc1 | tests first (2): E12b, open check, E13 mechanism, G6 (depth 2), C1 script, draft wording |
| 3893b19a | instrument repairs after R1 (G5's missing argument; E9's table reference; guarded follow-ups) |
| 0962cffb | the fix, ledger row `VAL-<train4-3>`, R1 records |
| f58cd9e2 | correction after R2: the replicate column's type as base left it when other runs move |
| 4e45c01e | merge of staging dd769671 (ledger to VAL-211; KIT-222's row after it, I06's corrections) |
| 92208fe5 | tests after the implementation audit: G1 population filter (I01 edit 3), E12a, E12d (I02), G7 |
| b8d9e56c | audit repairs under R/: base's type rule (I01 edit 1), no integer coercion in proteomics (I01 edit 2), I07 |
| 8414c68e | records after the audit (outcome.md, wording.md, base-tests.R) |

## Step-0 survey (before R1)

Tests that call the edits, metadata assignment, saves, imports or the drivers with a design holding `tech_reps`:
`helper-gui-entry.R:439-462`, `helper-module-ci-trx-tabs.R:86-98`, `helper-workflow-omics-residency.R` (manual,
corpus numbers, design order); `test-omics-design-core.R:400-403` (quarantined, Q-35), `:483` (lowest, in order);
`test-kit128-design-replicates.R` (lowest, in order); `test-prot-04b-design-server-shared.R:274-280` (lowest);
`test-workflow-gui-memory-lifetimes.R:296` (lowest); `test-lipid-13-design-builder-seams.R` (file skipped, Q-31).
No D1 or D2 within a key, no renumbering expectation, no metadata re-assignment of a marked run: no hit for the owner.
R2 confirmed: every one of these files `same`.

## Observed

### R1 (run 45507613511, e9309b9d on base): the defect in R

Every E row failed by assertion at its base value (59 failures; 3 errors, see instrument repairs). Observed base values
(lines `KIT222-OBS` in `r1/test-kit222.log`):

| Row | Ticket's base value | R on base |
|---|---|---|
| E1 four runs, 'first' | one key Ctrl_3, four runs | Ctrl_3 = {A_1, B_1, A_2, B_2} |
| E2 sorted, manual B=2, C=3, A=1 | Ctrl_2 = B + C | Ctrl_2 = {B_1, B_2, C_1, C_2} |
| E3 100-sample corpus, manual | 91 keys, nine named pairs | 91 keys: RPMI_17 = 17+18, RPMI_24 = 25+26, RPMI_42 = 42+43, RPMI_43 = 44+45, RPMI_46 = 47+48, Sera_11 = 11+12, Sera_13 = 14+15, Sera_28 = 30+31, Sera_46 = 46+47 (exactly the note's) |
| E3b trx summing | 91 units, nine sum four runs | 91 units; each of the nine sums four runs of two samples (e.g. RPMI_46: 116+131+145+196) |
| E3c 24 runs (G) | 12 keys | 12 keys (passes on both) |
| E4 lowest, three interleaved | B and C share 3 | Ctrl_3 = {B_1, B_2, C_1, C_2} |
| E5 manual 1 twice | ok TRUE, one key | ok TRUE, Ctrl_1 = four runs |
| E6 partial re-marking | tech_reps 1, 2, 2 | ok TRUE, tech_reps 1, 2, 2 |
| E7 lipid, metab, trx, proteomics builders, proteomics applet | one key | each: Ctrl_4 = {S2, S3, S4, S5} (A and B merged), applet's written design_matrix.tab too |
| E7 proteomics, number held | (head: shared text) | second edit accepted with the success text |
| E8 shared save on D1 | a payload | a payload (list) |
| E9 proteomics save D1, D2, OD-5 | accepted | R1: error in the test's table argument (instrument; repaired in 3893b19a); R2 base side: accepted, D1 Ctrl_3 = four runs, D2 Ctrl_3 = {A_1, A_2, C}, OD-5 Ctrl_1 = {A_1, A_3} |
| E9 proteomics checkpoint D1, D2 | accepted | error (no owner in the test; the assertion of the refusal class failed) |
| E10 lipid import D1, D2 | accepted | accepted |
| E11 resolver, RUV-III, lipid and metab fits, integration column, trx sum, peptide imputation | no condition | each returned (list, MArrayLM, character, PeptideQuantitativeData); none refused |
| E14 | one key, tech_reps 1, 1, 2, 2 | Ctrl_2 = four runs, tech_reps 1, 1, 2, 2 |
| E15 | one key | Ctrl_2 = four runs |
| E16 metadata route | payload, Trt_1 three runs | payload, Trt_1 = {A_3, C_1, C_2}; A_3 kept tech_reps 3 |
| E16b proteomics metadata route | payload | saved, Trt_1 = {S3, S4, S5} |
| consumers (builder's design) | merged | trx sum Ctrl_4 = S2+S3+S4+S5 = 14; DA block (lipid fit) and RUV-III sets: S2 to S5 one block; integration: 5 averaged columns, 4 of 6 names shared |
| G1, G4 | pass | pass |

Late rows, base values from R2's base side (3893b19a): E12b no message; open check absent (error); E13 nothing
recorded; G6 first violation at depth 2: lowest(r1, r2) then lowest(r1, r3).

### R2 (run 45509856651, 0962cffb against 3893b19a)

- 71 regression files (ticket section 12 grep, the GUI-entry files, the consumers' files): every one `same`.
- R CMD check: identical findings on both sides, 0 only on head.
- New file on head: 445 passed, 2 failed: G1 and G5, both because the replicate column was integer on head where
  base left it double (labels identical).
- C1: 54,024 G1 sequences; labels differ from base only in the 1,640 sequences where base merged (1,440 at four samples, 200 at fifty;
  all 'manual'); G6 depth 3 on head: no violation in 140,910 edits; base violates at depth 2.

### Correction pass (run 45519323914, f58cd9e2 against 3893b19a)

- The edit now keeps the column double when other runs' numbers move. 19 files that drive the edit, and check: all
  `same`, check identical. G5 passes. C1 G6 head: none.
- New file: 446 passed, 1 failed: G1, in 603 sequences of the test file (AUDIT-KIT-222-impl I01): 536 type-only
  (8 at two samples, 528 at three) and 67 random 50-sample 'manual' orders in which base merged, outside G1's declared
  population (the file's G1 had no base merge-free filter, unlike C1). In C1, 1,896 'manual' sequences (typed numbers in increasing order, out-of-row-order
  rows, e.g. rows A_1, B_1, A_2, B_2, B typed 1 then A typed 2) give identical labels but the column is integer on base
  and double on head: base moved no run there (nothing above the selection's highest number), while the head's ranks
  moved one. Type only; keys, partitions and every consumer's input values are identical.
- The exact repair is one condition: make the column double only when base would have shifted a run (any other run of
  the group numbered above the selection's highest number before the edit). The proteomics wrapper's integer coercion
  also moved the type against base, untested. Both were applied after the implementation audit (below).

### Audit-repair run (run 45574235919, 8414c68e against staging dd769671)

Authority: AUDIT-KIT-222-impl-20261125 and the coordinator's ruling (repairs after an independent audit, one run,
beyond the ticket's one correction; the owner told before it ran).
- 40 files (R2c's scope, `test-prot-design-log2.R`, the 18 integration files of I08 and the two integration-input files):
  every one `same`. R CMD check: identical findings, 0 only on head.
- The new file on head: 456 passed, 0 failed: G1 in full, G5, G7 (proteomics edit leaves base's type: double after a
  shift, integer otherwise), E12a and E12d (the key named on adopting a native archive and on opening a peptide
  checkpoint).
- The new file on base (`base-tests.R`, a worktree of dd769671 in the scripts job): every E row fails, E12a and E12d
  included (`r3a/base-tests.tsv`); G rows pass.
- C1: head merged 0; G6 head: no violation in 140,910 edits at depth 3. 52,336 of 54,024 sequences identical: the 1,640
  in which base merged, and 48 four-sample `manual` sequences (e.g. rows D_1, A_1, D_2, B_1, C_1, C_2, B_2, A_2 typed
  A=1, B=2, C=3, D=4) with labels identical to base and the replicate column of another type. INFERRED: base's shift
  condition is evaluated on base's intermediate numbers, which differ from the head's ranks during a multi-step
  sequence, so a per-step condition cannot reproduce base's type on every path. Not predicted (the audit: 52,384).
- F5 (I09) not done: the builder's selection observer (`R/mod_omic_design_builder.R:524-541`) re-sends each kept
  selection, `tech_rep_samples` included, after every design change, so a one-line clear after Assign would be
  overwritten in the same flush. A fix needs that observer changed as well: more than the register's one line.

## Acceptance criteria

| # | Status | Evidence |
|---|---|---|
| A1 | met: every E row failed on base by assertion at its base value; G rows passed (G5 errored on an instrument defect, repaired) | R1, R2 base side |
| A2 | met at the audit-repair run: every E row written passes on head (E12a, E12d included); E12c has no route (I10); E13 at the mechanism only. Earlier, at R2c: every written E row passed on the fix, but E12a and E12d were not written, so A2 was not met (it was recorded as met; AUDIT-KIT-222-impl I02); E13 is tested at the mechanism only. See the audit-repair run below | R2, R2c |
| A3 | met with the stated limit (coordinator's ruling, 2026-11-25, against stop rule 5; authority: AUDIT-KIT-222-impl-20261125, closure check I12): no head merge (G6 head none at depth 3; C1 head merged 0); every base merge absent on head; G1 in the test file passes in full; C1: labels and values identical in every base merge-free sequence; in 48 of 40,320 four-sample `manual` orders the replicate column's storage type differs (integer against double, values identical). Mechanism (hand trace of sequence 43986: rows D_1, A_1, D_2, B_1, C_1, C_2, B_2, A_2, typed A=1, B=2, C=3, D=4): the head's ranks compact the unmarked D_2 that base leaves at 3, so when D is marked, C (3) is above D's highest number (2) on head and the column becomes double, while on base D's runs are at 1 and 3, nothing is above, and it stays integer. The type reaches no number, file, key, label or record text; it reaches R-object digests (the lipid/metab saved DA fit identity, the proteomics base fit, the DIA product design id), so a fit saved by an earlier build is recomputed, with the same numbers, for those designs; no stored reference moves. No further run | r3a/c1 |
| A4 | met: every regression file `same`; check identical | R2, R2c, audit-repair run |
| A5 | not run (full suite after the audit, by the coordinator's order) | — |
| A6 | met except NEWS by decision: no NEWS (OD-6); ledger row naming IMP-036; `man/`, `NAMESPACE`, `DESCRIPTION` unchanged | diff |

## Deviations and inferences

- Tests were committed in three tests-only commits (e9309b9d, bcdd8bc1, 3893b19a) before the fix; R1 ran on the first,
  the late rows' base values come from R2's base side (3893b19a).
- E12a and E12d were not written at R2 on the stated reason that no builder at the branch point gives a native archive
  without a package change; that reason was wrong (AUDIT-KIT-222-impl I02): test fixtures at
  `test-workflow-native-session.R:5-22`, `:842-880`, `:965-1010` build both archives in code. They are written after
  the audit (92208fe5). E12c's premise is false (I10): no transcriptomics project can be reopened (`session_load` is
  configured for lipid and metab only, `R/utils_omic_gui_config.R:2477`, `:2695`; `workflowOmicsSessionClass` is
  lipid/metab, `R/utils_workflow_dia_sessions.R:1426`; no transcriptomics `owned_generation`), so OD-3's open message
  has no transcriptomics entry and nothing users have is left unproven there.
- E13 is tested at the mechanism only (the replay option records instead of refusing; the statement text), not by a
  full record replay (I03). The smallest full test, for later: one lipid record replay from an existing replay fixture
  with D1 rows in the design and `tech_rep_group` declared, expecting completion, `design_statement` naming the key and
  the warning class.
- The ticket's prose (sections 6, 7, 14 and G1's `manual` clause) said every in-order `manual` design is unchanged; R
  shows base merging in 1,440 of 40,320 four-sample row orders and 200 of 200 random 50-sample orders, as the plan
  audit's model had computed (Re-check 1, R1-03: 38,880 of 40,320 base merge-free). Example: rows D, A, B, C, typed
  A=1, B=2, C=3, D=4: B and C share 2. Corrected in ticket r1.3.
- The proteomics checkpoint refuses after the builder commit wrote its files (the builder's save refuses earlier); on
  the proteomics import path the refusal comes after the import's sidecars are read (INFERRED from the code).
- G6 in the test file is depth 2; depth 3 is in C1.

## Core-hours

R1 0.7, R2 28.6, correction 12.2, audit-repair run 20.0: 61.5 of the expected 107.7 (ceiling 501.33); R3 not run.

## Not checked

- Whether every class's design keeps `tech_reps` through every step (if a step drops it, the consumer check sees
  nothing there).
- A full replay of a record with a merged design; the four open sites end to end (E12a, E12d).
- Transcriptomics sessions on open (no entry found that adopts them through the named function).
- The full suite (R3).

## For the change report

- **What was wrong.** Marking technical replicates in any omic's design builder (and the proteomics workbook applet)
  renumbered the group's other runs, including runs already marked as another biological sample, so two biological
  samples could get one `tech_rep_group`. Every mode merged; Assign Metadata kept stale marks; proteomics accepted an
  unmarked run inside a marked sample's key. Shown in R (R1): the four-run example gives one key for two samples; the
  100-sample transcriptomics corpus gives 91 samples with the nine named pairs. Legacy since 33350d50 (2025-07-09).
- **Who.** Any user who marked technical replicates with 'Specify replicate number' or 'Use first selected sample's
  replicate', with 'Use lowest replicate number' when a sample had three or more runs not next to each other, with modes
  mixed, or who re-assigned marked runs; in proteomics also single-run samples. Users of the released `main`
  (proteomics DA block since 2026-03-14, peptide imputation key since 2025-12-11) are not notified (OD-8).
- **Code.** Base: `R/func_omics_design_shared.R:629-655` (the shift `:643-648`), `:610-623`, `:724-747`;
  `R/utils_workflow_dia_state.R:1229-1260`, `:1416-1424`; `R/utils_workflow_omics_design.R:255-290`. Fix (f58cd9e2):
  `applyOmicsDesignTechRepAssignment`, `findOmicsDesignReplicateKeyCollisions`, `refuseOmicsDesignReplicateKeyCollision`,
  `omicsDesignSaveRefusal`, `noteOmicsDesignReplicateKeyCollisionOnOpen` (`R/func_omics_design_shared.R`); the check in
  `resolveProtTechnicalReplicateBlock` (`R/func_prot_norm_optimization_helpers.R`); the four open sites
  (`R/utils_workflow_dia_sessions.R`); replay (`R/utils_workflow_analysis_record.R`).
- **What changed in results.** Designs with a merged key: now refused at save, import and at every step that computes
  numbers, or (by construction) not produced by the builder; U8's corpus gives 100 samples instead of 91. Merge-free
  designs: identical labels and partition to base on every G1 sequence (C1); the replicate column's type differs in
  48 four-sample manual sequences after the audit repairs (open item; was 1,896 before them). IMP-036's stated limit (2026-11-17) understated the defect and is superseded.
- **Evidence.** R1 45507613511; R2 45509856651; correction 45519323914; C1 in `r2/c1`, `r2c/c1`. The ticket's and the
  audit's Python models predicted every base value R showed, including the 1,440 of 40,320 four-sample `manual` orders.

## Past results

How a user checks a saved design (in the app's messages, OD-6): in `Study_report/design_matrix.tab`, count the distinct
`tech_rep_group` values in each group; fewer than the biological samples put in that group means two samples were
analysed as one. A transcriptomics analysis can be checked only from its saved pre-sum design file. KIT-209's U8 and
D3 figures analysed 91 units (9 merged pairs).
