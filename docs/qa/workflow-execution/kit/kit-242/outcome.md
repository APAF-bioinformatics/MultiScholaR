# KIT-242: outcome (branch `kit-242`)

Ticket r1.4. Branch from `v0.9-staging` e0747fc198817596e3e66397a87e06b2f869c6a7. G0 and C4b are recorded on branch
`kit-242-g0` (records at 7ea50bc1, `docs/qa/workflow-execution/kit/kit-242/`): 2.76 core-hours.

## Tests-only commit 17488ccf: the option-free rows (ticket section 7, r1.4)

`tests/testthat/test-kit242-normalisation-fallback.R`, nothing under `R/`. Rows, lipid and metab: E2 (quantile on a
sample with one value stops, `omics_normalisation_failed`, `normaliser_error`, limma's error as `parent`; and through
the declaration, no normalised state committed), E3 (an assay with no sample columns stops, `malformed_assay`), E16
(legacy RUV-III-C with k = -1 for A2 stops, `omics_ruv_failed`, `invalid_k`), E10 (the record's v2 DA entries carry
`fdr_method`, "BH" with `fdr_fallback_reason` where the manifest records a q-value fallback), NA1 to NA3 (the record
states between-sample normalisation, ITSD and RUV per assay, FALSE for an assay the step left unchanged; every assay
the RUV record names is in the RUV state). Seeds 242000 + row number. None runs the correlation filter, so they
assume nothing about KIT-255's defect.

## R1 (run 45602184961, 9.9 job-minutes, 0.66 core-hours): observed

Every one of the 14 tests fails on base (10 expectations passed, 46 failed, 2 errors), at the base values G0 observed:

| Row | Base value observed |
|---|---|
| E2 | no error is raised (G0 C2: a warning and the assay unnormalised); the error expectation fails, then the fields are NULL; `conditionMessage(NULL)` errors, which ends the test before its declaration arm (that arm runs only once the error exists) |
| E3 | no error is raised; fields NULL |
| E16 | no error is raised (G0 C12: the assay is dropped with a logger line); fields NULL |
| E10 | `fdr_method` is NULL for both entries; in this draw A1 fell back too (as in G0 C7) |
| NA1, NA2 | the per-assay outcome is absent for A1 and A2 |
| NA3 | the per-assay outcome is absent, and A2, named by the record, is absent from the RUV state |

Stop rule 1 (a row passing on base) did not fire.

## Not checked

- These rows at head: they wait for the fix, which waits for the plan auditor's re-check of r1.4.
- P0: it runs on a scratch commit carrying the fix (ticket 8.2), so it cannot run before the fix.

## Part A (ticket r1.5): R1b-A, fix, R2a

Commits on `kit-242`: 4da0750c (merge of staging 7327607f), 17d1ecc4 (tests-only R1b-A: Part A's rows; E16 and NA3 move
to Part B), dc75de69 (Part A's fix and NEWS).

### R1b-A (run 45608276744, 13.9 job-minutes, 0.93 core-hours): observed

Every E and NA row fails at the base value G0 observed (85 expectations passed, 113 failed, 5 errors): E1 no
`omics_normalisation_sample_excluded` and A2 unnormalised; E2, E3 no error (the 2 errors per class are `err$assay` on the
returned object); E4 to E7 the per-assay fields absent; E8 no log line or notice naming A2 and B_6; E9 "- ITSD applied:
Yes"; E10 no `fdr_method`; E11 no marker; E13 the guard's `omics_cyclic_loess_too_few_features`; E17 scale all missing,
cyclic loess keeping B_6, quantile limma's error; E26 B_6 kept; NA1, NA2 the per-assay outcome absent. The guards G1,
G2, G3, G7 pass. Stop rule 1 did not fire.

### R2a (run 45610857369, base 17d1ecc4, head dc75de69, 138 files, 9 shards per side, check; 746.5 job-minutes, 49.77 core-hours)

Observed:
- `test-kit242-normalisation-fallback.R`: fixed, 247 passed, 0 failed at head.
- Check: findings identical on base and head (6; none only on head).
- 123 files `same`; 11 `new-failure`; 3 `message-changed`.

| File | Verdict | What moved | Class |
|---|---|---|---|
| `test-gui-entry-lipidomics.R` | new-failure (35/0 to 25/10) | lane captures after the ITSD step (log, record, files) | predicted (8.1, OD-2 (b)); by-file question (ii) OPEN |
| `test-gui-entry-metabolomics.R` | new-failure (105/0 to 75/30) | metab_lc, metab_gc and metab_combined captures | predicted class (8.1 named metab_lc only); (ii) OPEN |
| `test-gui-entry-proteomics.R` | new-failure (189/4 to 174/19) | prot_dia_limpa: the normalisation text "State saved as" spacing | **defect of the fix** (an extra `paste()` element in `computeDIAReferenceProteinNormalisation`'s text changes spacing when nothing is removed); whether the provenance marker also moves this lane is to be re-checked after the correction |
| `test-kit172-downstream.R` | new-failure (84/0 to 70/14) | T10: one more warning and notification (the new ITSD notice on KIT-172's fixtures, where ITSD matches nothing, K242-P01) | predicted class (8.1, KIT-172 fixtures) |
| `test-kit172-reopen.R` | new-failure (14/0 to 12/2) | T9: a project saved by the base commit replays with mismatches in the ITSD, log2 and normalised states' products and summaries | **unpredicted**: the ITSD step's state arguments and product now carry per-assay fields, and replay compares products and summaries, so an older record no longer replays exactly; N05's "no marker on the ITSD step" does not hold |
| `test-metab-02w-…-characterization.R` | new-failure (16/0 to 13/3) | the all-missing-sample expectation: no warning; the fixture's only assay has Sample_2 empty, so Sample_2 leaves the design (E26's rule) | predicted (8.1 (i)); by-file question (i) OPEN |
| `test-omic-norm-module.R` | new-failure (639/0 to 583/56) | provenance digests (the marker in the step's `fixed`) | predicted (8.1) |
| `test-workflow-lipid-declarations.R`, `test-workflow-metab-declarations.R` | new-failure (259/0 to 257/2 each) | `itsd_applied` TRUE to FALSE on a fixture whose ITSD matches nothing | predicted (8.1) |
| `test-workflow-omics-worker.R` | new-failure (289/0 to 275/14) | the raw declarations fixture (duplicate feature ids, duplicate resolution not run) now stops at between-sample normalisation: "Error converting tibble to matrix: duplicate 'row.names' are not allowed" | **unpredicted**: 2.1 called the matrix-conversion row unreachable; a duplicate id reaches it |
| `test-workflow-operation-registry.R` | new-failure (1095/2 to 1094/3) | KIT-025's scan for random-number calls matches the text "sample(s)" (`\bsample\s*\(`) in `computeDIAReferenceProteinNormalisation` | **defect of the fix** (wording) |
| `test-prot-norm-core.R`, `test-workflow-omics-binding-steps.R`, `test-workflow-omics-declarations.R` | message-changed (failing on base already) | the frozen KIT-073/KIT-094 references: the new marker parameter rows | predicted (8.1); by-file question (iii) OPEN |

Stop rule 2 fired: stored references the owner has not approved by file moved (expected, the by-file questions being
open), and three unpredicted verdicts. No correction was made.

Inferred (not executed): the correction would (1) reword "sample(s)" and keep the text's spacing unchanged when nothing
is removed; (2) put the ITSD step's new per-assay arguments and product fields under a replay marker, so an older record
replays with the earlier ones (a change to N05, for the auditor); (3) decide whether a duplicate-id assay at
between-sample normalisation is "malformed" (stop, moving the worker test's expectation by file) or keeps today's skip.

### A3 (Part A's methods)

`git grep -n "warning(\|log_warn(" -- R/func_lipid_s4_normalization_methods.R R/func_metab_s4_transformation_methods.R`
at dc75de69: in `normaliseBetweenSamples`, every warning followed by a return of an unnormalised assay is in the
earlier-rule branch (`current_rule` FALSE: a replay of an older record) or inside `malformed()`, which stops under the
current rule; in `normaliseUntransformedData`, the remaining skips are OD-2 (b)'s decided handling (`notApplied`,
recorded and noticed). Not checked: the RUV methods (Part B).

## Not checked (Part A)

- A separately started worker's log; the tab's notice on the proteomics routes (Part C's samples are recorded in the
  limpa route's step outcomes and text; the `normalization-ruv` route's warning is not captured into its record).
- The study-parameters writer's per-assay ITSD lines (in the file list; not written in Part A).
- The note OD-4 (c) asks for wherever an older record is read (not written in Part A).

## Part A's correction pass (the coordinator's rulings of 2026-11-26; one pass, no second)

Rulings: item 4, the protein normalisation text byte-identical to base when nothing is removed; item 3, "samples"
instead of "sample(s)" (the KIT-025 scan untouched); item 1, the ITSD step's per-assay record under a replay marker
(`itsd_outcome`), correcting N05; item 2, an assay with duplicate feature ids at between-sample normalisation is a
malformed assay and stops, naming the cause and the remedy (draft `wording.md`).

Commits: a5a87d7a (tests-only: E3b, E11b, the text pin, E8's wording, and the declared edit of
`tests/testthat/test-workflow-omics-worker.R:15-23`, which pinned the earlier skip on the raw fixture's duplicate ids).

### Predictions, written before any run of this pass

**R1c** (the new rows, on a5a87d7a = dc75de69's fix without the correction; files `test-kit242-normalisation-fallback.R`
and `test-workflow-omics-worker.R`): fail: E8 (4 tests, the wording "samples"), E3b (2: the refusal's message does not
name duplicated identifiers or 'Resolve Duplicates'), E11b (2: no `itsd_outcome` marker), the text pin (1: no
`diaProteinNormalisationText`); the worker test's between-sample check (2 expectations: the message); every other
test of both files passes as in R2a's head.

**R2a-2** (R2a's scope again, base 17d1ecc4, head = the correction commit, 138 files, 9 shards per side, check): the
four unpredicted verdicts of R2a gone: `test-kit172-reopen.R` same; `test-workflow-omics-worker.R` same (base passes
the old loop, head the edited loop); `test-workflow-operation-registry.R` same; `test-gui-entry-proteomics.R`
prot_dia_limpa back to same (if the provenance marker moves no capture of that lane; if it does, that is a stored
reference of 8.1's KIT-073/094 class, listed, not a new verdict). `test-kit242-normalisation-fallback.R` fixed (every
row passes at head). Everything else as R2a: the predicted moves (`test-gui-entry-lipidomics.R`,
`test-gui-entry-metabolomics.R` with metab_lc, metab_gc and metab_combined, `test-kit172-downstream.R` T10,
`test-metab-02w-…`, `test-omic-norm-module.R`, the two declarations files, and the message-changed
`test-prot-norm-core.R`, `test-workflow-omics-binding-steps.R`, `test-workflow-omics-declarations.R`), the rest same;
check identical. A verdict outside this: stop and report, no second correction.

### R1c (run 45621160740, a5a87d7a; 21.6 job-minutes, 1.44 core-hours): observed

As predicted: E8 (4 tests, the wording), E3b (2), E11b (2: no `itsd_outcome` marker; the earlier and current journeys'
fields), the text pin (1: no `diaProteinNormalisationText`) and the worker test's between-sample check (2 expectations:
the message) fail; every other test of both files passes (267 + 275 passed). Stop rule 1 did not fire.

Correction commit: 25fed307.

### R2a-2 (run 45625379503, base 17d1ecc4, head 25fed307, 138 files, 9 shards per side, check; 815.0 job-minutes, 54.33 core-hours): observed against the prediction

- **As predicted:** `test-kit172-reopen.R` same (T9: the project saved by the base commit replays exactly);
  `test-workflow-operation-registry.R` same; `test-gui-entry-proteomics.R` same (prot_dia_limpa); the new file fixed
  (290 passed, 0 failed); check identical (6 findings on both, none only on head); the predicted moves as in R2a for
  `test-gui-entry-lipidomics.R` (10), `test-gui-entry-metabolomics.R` (30), `test-kit172-downstream.R` (14),
  `test-metab-02w-…` (3), `test-omic-norm-module.R` (56), and the message-changed `test-prot-norm-core.R`,
  `test-workflow-omics-binding-steps.R`, `test-workflow-omics-declarations.R`; 126 files same.
- **Outside the written prediction (stop, as ruled; no second correction):**
  1. `test-workflow-omics-worker.R` is `count-changed` (passed 289 to 277, 0 failures), not `same`: every expectation
     passes; the declared edit skips the six later checks of the refused step per omic. A wrong verdict label in my
     prediction, not a failure.
  2. `test-workflow-lipid-declarations.R` and `test-workflow-metab-declarations.R` each move one more expectation than
     in R2a (3 against 2): KIT-031's and KIT-032's registration tests pin the ITSD declaration's `fixed` as
     `list(method = "ITSD")`, and the ruled marker adds `itsd_outcome`. A stored expectation of the marker class
     (8.1's provenance row) that my prediction did not list.

### The stored expectations that moved, by file and key (to the owner with the three by-file questions)

Full list, one row per moved expectation: `r2a2-45625379503/moved-expectations.tsv` (file, test, expectation, key).
No stored value was edited.

| File | Moved | Keys |
|---|---|---|
| `test-gui-entry-lipidomics.R` | 10 | lane lipid_canonical, steps da, da_load, session_export, summary_copy, summary_export, summary_refusal, summary_report, summary_save, and the norm log output (2) |
| `test-gui-entry-metabolomics.R` | 30 | lanes metab_lc, **metab_gc**, **metab_combined** (the last two not named in 8.1), each at steps da, da_load, session_export, summary_copy, summary_export, summary_refusal, summary_report, summary_save, and the norm log output (6) |
| `test-kit172-downstream.R` | 14 | T10 lipid and metab: the first warning's text and the warning count (in process, managed); the count of warning notifications |
| `test-metab-02w-normalise-between-samples-characterization.R` | 3 | "preserves all-NA sample fallback": the warning, `design_matrix$Run`, the assay table (Sample_2 has no value in the only assay, so it leaves the design) |
| `test-omic-norm-module.R` | 56 | (b) lipid and metab, RUV skip and ITSD off (5 each); RUV manual with a manual ITSD selection (12 each); (b, c) RUV automatic (11 each): provenance digests |
| `test-workflow-lipid-declarations.R` | 3 | the ITSD declaration's `fixed`; the journey's `itsd_applied` outcome; the report overlay's `itsd_applied` |
| `test-workflow-metab-declarations.R` | 3 | the same three for metab |
| `test-prot-norm-core.R` | (message-changed; base already failing) | KIT-111 W1 (a) frozen references: the new marker parameter rows |
| `test-workflow-omics-binding-steps.R` | (message-changed; base already failing) | KIT-094 snapshot: the steps' parameter lists |
| `test-workflow-omics-declarations.R` | (message-changed; base already failing) | KIT-073 frozen dump |
| `test-workflow-omics-worker.R` | declared edit (`:15-23`) | the between-sample step on the raw fixture refuses (duplicate ids) |

Budget: 210.78 expected / 782.0 ceiling core-hours, as the coordinator told the owner. Used to date: G0 2.76, R1a 0.66,
R1b-A 0.93, R2a 49.77, R1c 1.44, R2a-2 54.33 = 109.89.

## Part B (the coordinator's rulings of 2026-11-26; base 0020144c)

Ceiling: 901.33 core-hours for the ticket, accepted by the owner (ledger VAL-329 on origin/v0.9-staging 176f9d3f).
Rows fixed at the branch point: R1b-B 2.0 expected; P0 28.0 expected; R2 59.33 / 120.0; R3 77.7 / 186.67. Order: tests
first, R1b-B, the fix, P0, then the moved-expectation list and the three by-file questions to the coordinator before R2.

Tests-only commit fb5eee2b (nothing under `R/`): Part B's rows at the end of `test-kit242-normalisation-fallback.R`
(E14, E15, E16, E18, E19, E20, E22, E23, E24 with the tab managed and in process, E25, E25e (VAL-290) with the tab in
process, NA3, E11c, the row of ruling (7), G4, G5, G8, G9; seeds 242000 + k with k = 14, 15, 16, 18, 19, 20, 22, 23, 24,
25, 125, 127, 12, 44, 48, 45, 49), and the two ITSD `fixed` pins re-pinned (VAL-303 (2)). `kit242Tab` also returns
`ruv_complete` and the tab's RUV results.

### Predictions for R1b-B, written before the run (fb5eee2b: base code 0020144c plus the tests; the new file only)

- Every row written before Part B (Part A, the correction pass) passes, as on R2a-2's head.
- Each of these Part B tests fails, per class (17 per class, 34 in all), at the base value named: E14 (A2 absent from
  the returned object; no `omics_ruv_not_applied`), E15 (no classed warning), E16 (no error), E18 (the step completes
  with A2 absent), E19 (`omics_ruv_insufficient_controls` stops the step), E20 (A2 has no complete feature, so the typed
  k is refused and stops: `ruv_k_above_limit`; the ticket's "refused per assay" is the automatic rule's form of the same
  refusal), E22 and E23 (A2 absent), E24 (A2 absent; no per-assay record), E24 tab managed and in process (no "RUV-III
  was not applied to assay A2" line), E25 (no error: A2 dropped), E25e (both assays dropped and the input returned;
  `ruv_applied` TRUE), E25e tab (the log says "RUV-III correction applied"), ruling (7) (`ruv_complete` TRUE), NA3 (A2
  absent from the RUV state), E11c (no `ruv_fallback` marker).
- G4, G5, G8 and G9 pass (8 tests).
- Stop rule 1: a Part B E, NA, E11c or ruling (7) row passing on base, or a guard failing, stops before the fix.

### R1b-B (run 45714743051, fb5eee2b; 23.3 job-minutes, 1.55 core-hours): observed

- As predicted: each of the 17 Part B tests per class fails at its base value (82 failed expectations, 4 errors; the
  errors are `expect_match` on the base's absent `ruv_not_applied_text` in E24 and E25e). E18, E22, E23: A2 absent; E19
  and E20: the step stops; E22 records A2 `dropped`; E15: no classed warning.
- G4, G5, G8, G9 and every row written before Part B pass (370 passed). Stop rule 1 not triggered. Tables in
  `r1b-b-45714743051/`.

### Fix 5c856dff; prediction for R1d (the new file on the fix, before P0), written before the run

The smallest falsifying test of the fix before P0's 186 files: `test-kit242-normalisation-fallback.R` on 5c856dff (the
tests-only commit plus the fix), 1 job, ceiling 60 job-minutes. Prediction: every test passes (Part A's rows, the
correction pass's rows, Part B's 42 tests including G4, G5, G8, G9).

### R1d (run 45721168959, 5c856dff; 13.9 job-minutes, 0.93 core-hours): observed

As predicted: 512 passed, 0 failed, 0 skipped. (Dispatcher run 45721099733 only launched it.)

### P0, declared before the run

Head-only job over the 186 files of `r2-files.txt`, 12 shards, no check, on scratch commit `kit-242-p0` = 0020144c plus
the fix 5c856dff (no test edits: the tests-only commit fb5eee2b is not in it, so the new file runs in its Part A form
and the two ITSD `fixed` pins at their old value). Expected 420 job-minutes (28.0 core-hours), ceiling 840 (56.0).
Purpose: list every stored expectation that fails on the fix, then classify each by the marker or rule that moves it
and whether it already failed on runners (R2a-2's tables for Part A's moves and the files' base failures; the staging
train runs where they cover a file; otherwise stated as not known before R2). No prediction of the list is made (P0 is
the measurement); an instrument failure (a shard that does not complete) is reported, not counted.

### P0 (run 45723354055, 0c8778d9; 539.2 job-minutes, 35.95 core-hours; ceiling 840): observed

186 files: 152 ok, 34 failed, none timed out or crashed. Full list, one row per failing expectation with its class:
`p0-45723354055/moved.tsv`. Classes (observed by matching each failure to the runner tables named):

| Class | Expectations | Files |
|---|---|---|
| Moved by Part B (new) | 16 | `test-workflow-omics-worker.R:25-31` (7 each for lipid_ruv, metab_ruv): the RUV step on the raw fixture now stops; `test-module-ci-metab-norm.R:202`, `:204` (MCI-015.3): a failed Manual search now stops |
| Defect in the fix (Part B) | 1 | `test-workflow-operation-registry.R:257` (KIT-025 scan): the text "every sample (%s)" in `omicsRuvSearchOutcomes` matches `sample\s*\(` |
| Moved by Part A (R2a-2's list) | 79 | `test-omic-norm-module.R` 56 (call sites `:1560`-`:1585`), `test-kit172-downstream.R` 14, `test-metab-02w-…:226-240` 3, the two declarations files 3 each (`:47`, `:138`, `:188`) |
| Already failing on staging train heads | 951 | proteomics, enrichment, summary, DA modules, KIT-073/094 fixtures, GUI-entry lanes (lipid 16, metab 48, proteomics 101), others (`moved.tsv`) |
| Transcriptomics corpus-manifest digest | 22 | `test-kit129-run-path.R`, `test-kit224-…`, `test-module-ci-trx-norm.R`, `test-trx-da-without-ruv.R`, `test-trx-ruv-designs.R`, `test-workflow-trx-declarations.R`: fail on R2a-2's base 17d1ecc4 too, pass on trains 15/16 (environment, not a KIT-242 rule) |
| Runner library | 2 errors | `test-tech-reps-limma.R` (devtools), `test-workflow-checkpoint-harness.R` (lobstr) |

Inferred, not observed: inside expectations that already fail, `ruv_fallback` also moves the RUV step's `fixed` list
(KIT-073 dump, KIT-094 snapshot, the GUI-entry lanes' RUV steps) and the provenance of `test-omic-norm-module.R`'s manual
and automatic sessions; P0 cannot separate these (the expectation fails either way). The worker's message ("the control
features ... must be feature ids, positions or a logical vector") is inferred to come from the Manual search returning
no per-assay control set on the raw fixture (the optimiser then hands the whole list); the message names the wrong cause.

## Not checked (Part B so far)

- R2 and R3 (owner's answer on the stored values first); the declared design through the tab; KIT-188/KIT-262 journeys.
- Whether any expectation already failing on staging hides a numeric move by Part B (R2 base-against-head decides).

## For the change report (Part B)

Legacy lipid/metab RUV-III-C dropped an assay it could not correct (and returned the input when every assay failed)
while recording `ruv_applied` TRUE; on main since 2025-05-06 (metab) and 2026-01-18 (lipid). Now by cause (VAL-242,
VAL-290), with marker `ruv_fallback`; evidence R1b-B 45714743051, R1d 45721168959, P0 45723354055.

## After P0: the coordinator's rulings (2026-11-26), no run

**Cost against the lines.** R1d (45721168959) was not a planned row: run as the smallest test that could falsify the fix
before spending P0's 186 files on it (0.93 core-hours). P0 used 35.95 core-hours against its line of 28.0 expected and
56.0 ceiling. Used to date: 148.31 of 901.33 (VAL-329).

**Correction 1 (accepted).** The no-usable-controls notice reads "observed in every sample: <reason>" (was "every
sample (<reason>)", which KIT-025's scan matched). Commit 4c8e3d87.

**Correction 2: the cause, read (observed in the code, not run).** The worker fixture's assays each hold two duplicated
feature ids (`tests/testthat/helper-workflow-lipid-declarations.R:29-31`: `ids[c(3L, 7L)]` repeated; metab twin). The
control selection converts each assay with `tibble::column_to_rownames` (`R/func_lipid_s4_ruv_control_methods.R:148-161`;
metab `R/func_metab_s4_design_control_methods.R:235`), which fails on duplicated row names, so the assay is skipped and
dropped from the returned list (`:348`; metab `:431`). The Manual search then takes the whole list as that assay's
control set (`R/func_lipid_norm_ruv_optimization_helpers.R:421-425`; `R/func_metab_norm.R:944-948`, code unchanged by
KIT-242) and reports success. So the search really returns no per-assay control set: the stop is Q6.2 (c)'s (a failed
search), and only the message was untruthful (it came from `omicsRuvCorrectByCause`'s control-type check). Repaired in
4c8e3d87: `omicsRuvSearchOutcomes` stops such an assay first (reason `search_error`), with the message (draft for the
wording pass): "RUV-III was not run: the control-feature search found no control features for assay <A>, because it
could not read the assay. The assay has <n> duplicated feature identifiers (<first three>). Run 'Resolve Duplicates' in
the Quality Control tab, then run the normalisation again." (the second and third sentences only when duplicates are
present). Test row E25m (lipid and metab).

**Correction 3: the two declared edits, prepared, not committed (to the owner with (i) to (iv)).**

`tests/testthat/test-workflow-omics-worker.R`, after the between-sample branch (`:22`), insert:

    # KIT-242 Part B (VAL-242 Q6.2 (c); declared edit, owner <date>, ledger <row>): on the raw fixture the control
    # selection cannot read an assay with duplicated feature ids, so the RUV step stops, naming the cause and the
    # remedy; lines 25-31 pinned the earlier completion (the assay dropped by the correction).
    if (identical(step$operation, paste0(binding$prefix, "_ruv"))) {
      expect_identical(result$status, "failed", info = paste(omic, step$operation))
      expect_match(result$error %||% "", "duplicated feature identifiers.*Resolve Duplicates", info = omic)
      next
    }

(lines 25-31 then apply to the other steps only, unchanged.)

`tests/testthat/test-module-ci-metab-norm.R:200-204`, old:

    failed <- run_ruv("manual", function(...) failed_results, manual_percentage = 1)
    expect_true(is.na(failed$best_k$LCMS_Pos))
    expect_null(extractCtrlPerAssay(failed_results)$LCMS_Pos)
    expect_identical(failed$best_k$GCMS, 2L)

new:

    # KIT-242 Part B (VAL-242 Q6.2 (c); declared edit, owner <date>, ledger <row>): a failed Manual search stops the
    # step, naming the assay, where LCMS_Pos was dropped and GCMS corrected alone.
    failed <- run_ruv("manual", function(...) failed_results, manual_percentage = 1)
    expect_null(extractCtrlPerAssay(failed_results)$LCMS_Pos)
    expect_false(failed$state$result)
    expect_true("ERROR: RUV-III was not run: the search for assay LCMS_Pos failed: too few controls" %in% failed$state$logs)
    expect_false(failed$norm_data$ruv_complete)
    expect_null(failed$best_k)

Note for the owner: in the real Manual search `success = FALSE` arises only from an error caught at
`R/func_lipid_norm_ruv_optimization_helpers.R:451-466` (metab twin), so a failed Manual search is a failed computation;
the stub's text "too few controls" is the test's.

**Request environment (item 4).** P0's request set `env: {}` and `extra_packages: []`. The trains set
`MULTISCHOLAR_TRX_MANIFEST_NOTE=1` (trains 15 to 19, 22 q) and, where their files need them, `extra_packages`
`["limpa", "mixOmics", "sva"]` and `file_timeouts` of 120 for the GUI-entry and prot-norm-core files. The 22
transcriptomics corpus-manifest failures are therefore attributed to P0's request, not the code (inferred from the
requests; R2 confirms). R2 and R3 will set the trains' `env`, `extra_packages` and `file_timeouts`. Of the two
runner-library files: `lobstr` (`test-workflow-checkpoint-harness.R:187`) is not in `packaging/multischolar.lock.json`, so
no `extra_packages` entry can supply it; `test-tech-reps-limma.R:6` calls `devtools::load_all(".")` (devtools is in the
lock). No train request of 2026-11-26 lists either file. R2 and R3 leave both out, as the trains do, with this reason;
neither touches a KIT-242 path.

## The owner's answers (VAL-352, origin/v0.9-staging 91713498) and R2's predictions, written before R2

Tests-only commit 51b77bb2 (nothing under `R/`): `test-metab-02w-…` all-NA case re-pinned to OD-1 (c) (VAL-352 (1));
the declared edits of `test-workflow-omics-worker.R` (RUV branch) and `test-module-ci-metab-norm.R` MCI-015.3 (VAL-352
(5)). Answers (2) to (4): nothing edited.

**R2's shape (sent only on the coordinator's word, after the one merge of staging once train 25 has landed).** Base: every
Part B test change (fb5eee2b, the E25m row of 4c8e3d87, 51b77bb2) on the code before Part B's fix (staging merge plus
Part A); head: kit-242's head after the merge. Files: `r2-files.txt` recounted at the merge, without
`test-workflow-checkpoint-harness.R` (lobstr is not in the lock) and `test-tech-reps-limma.R` (`devtools::load_all`), as
the trains; `env` `MULTISCHOLAR_TRX_MANIFEST_NOTE=1`; `extra_packages` `limpa`, `mixOmics`, `sva`; `file_timeouts` 120
for the GUI-entry files and `test-prot-norm-core.R`; check on both sides. Rows fixed at the branch point: 59.33 / 120.0.

**Predictions by file (verdicts of `compare.tsv`).**

| File | Base | Head | Verdict |
|---|---|---|---|
| `test-kit242-normalisation-fallback.R` | the Part B tests fail (17 per class as R1b-B) and E25m fails by value (c93c4848: the function is absent; no error) | all pass | fixed |
| `test-workflow-omics-worker.R` | the new RUV branch fails for lipid_ruv and metab_ruv (the step completes on the code before the fix) | all pass | fixed |
| `test-module-ci-metab-norm.R` | MCI-015.3's edited block fails (the run completes with LCMS_Pos dropped) | all pass | fixed |
| `test-metab-02w-normalise-between-samples-characterization.R` | passes (Part A's code is on the base) | passes | same |
| `test-workflow-lipid-declarations.R`, `test-workflow-metab-declarations.R` | `:138`, `:188` fail (Part A, re-record list); `:47` passes (VAL-303 re-pin) | the same two, identical messages (ITSD outcomes; no RUV value in them) | same |
| `test-omic-norm-module.R` | the 56 fail (Part A) | see below | message-changed, exactly the keys below |
| `test-workflow-omics-declarations.R` (KIT-073) | fails as on staging (registry, values, states, correlations, DA) | `current$registry` (the `lipid_ruv` and `metab_ruv` declarations' `fixed` gain `ruv_fallback`) and `current[[omic]]$values` for lipidomics and metabolomics (the RUV step's value gains `ruv_applied`, `ruv_not_applied_reason`, `ruv_not_applied_text`, `ruv_samples_excluded`, `ruv_features_excluded`) change their actual values; states, correlations and DA identical | message-changed, those 3 expectations |
| `test-workflow-omics-binding-steps.R` (KIT-094) | the `report_sources` and `records` parts fail as on staging | the `records` part changes its actual value (the RUV step's record gains the parameter row `ruv_fallback` and the outcomes `ruv_applied_per_assay`, `ruv_not_applied_reason_per_assay`); `report_sources` identical; `registry` passes on both (it lists no `fixed`, test file `:57-63`) | message-changed, the `records` expectation |
| `test-gui-entry-lipidomics.R`, `test-gui-entry-metabolomics.R` | fail as on staging | identical: the lipid and metab lanes set `norm-ruv_mode` "skip" (`helper-gui-entry.R:615`), so no legacy RUV step, marker or value enters them | same |
| `test-workflow-operation-registry.R` | staging's 2 failures | the same 2 (the KIT-025 scan passes after correction 1) | same |
| `test-kit172-downstream.R` | 14 (Part A) | 14, identical (T10 is the between-sample `scale` warning; no RUV) | same |
| the transcriptomics files (`test-kit129-run-path.R`, `test-kit224-…`, `test-module-ci-trx-norm.R`, `test-trx-da-without-ruv.R`, `test-trx-ruv-designs.R`, `test-workflow-trx-declarations.R`) | with the trains' `env` the 22 manifest digests pass; any staging failures stay | the same | same |
| every other file | as on staging after the merge | the same | same |
| check | | | identical findings on both sides |

**`test-omic-norm-module.R`, by reading (the frozen captures `:1560`-`:1585`, compared at `:1549`).** The test stubs the
search and replaces `ruvIII_C_Varying` with a stub that returns the object with its assays unchanged
(`test-omic-norm-module.R:283-289`). On the head the kernel therefore records both assays `uncorrected`, `ruv_applied`
FALSE, and the tab's closing line becomes "RUV-III not applied: no assay could be corrected (see the warnings above);
the data continue without RUV" (`omicsRuvAppliedText`) in place of "RUV-III correction applied". The provenance
captures carry each step's `fixed` (`omicNormProvenance`, `:244-249`), so the RUV step's `ruv_fallback` enters them. The
states are the S4 states (unchanged by the stub), so they do not move; no notification is added (no notice text).
Expectations whose actual value changes from base to head (both omics):
- `log`: manual run, automatic run (4 keys);
- `provenance`: manual run, manual skip_correlation, manual reset, automatic run, automatic final_plot, automatic
  final_plot_failing (12 keys);
- `outputs`: manual run, manual skip_correlation, automatic run, automatic final_plot, automatic final_plot_failing (10
  keys; the `norm_log` output holds the cumulative log, `R/mod_omic_norm_panels.R:1136-1145`).
- Not settled by reading: manual `reset outputs` (2 keys): whether the reset clears the log before the capture. Not
  predicted; observed and recorded: changed means the log survives the reset, identical means it is cleared.
- Identical on both sides: the skip sessions (10), every `states` and `final_pca` key, and `normdata` (passing on both:
  no assay is kept, so the RUV results are not annotated).

Any other verdict (a new failure, a file fixed that is not listed, a count change in a file not listed) is outside the
prediction: stop and report.

**Change to E25m (c93c4848), on the base only.** On the code before Part B `omicsRuvSearchOutcomes` does not exist, so the
row now fails by value there (`expect_true(is.function(...))`) instead of erroring; what it asserts on the head is
unchanged.

## Merge of staging aa119308 (afb27738) and R2 restated, written before R2

**Merge (once, on the coordinator's word).** Two files conflicted, nothing else:
- `R/utils_workflow_analysis_record.R`, the marker list. Merged order: missing_value_comparison, enrichment_background,
  intensity_filter_samples, duplicate_resolution_samples, ruv_feature_levels, treat_moderation, ruv_control_grouping,
  ruv_k_bound, intensity_filter_missing, duplicate_resolution_rank, correlation_exclusion, normalisation_fallback,
  protein_normalisation_fallback, itsd_outcome, ruv_fallback.
- `R/utils_workflow_saved_exports.R`, the session list: staging (KIT-255) made the list a variable and adds
  `correlation_exclusions` when the state carries any; KIT-242 adds `itsd_outcome_per_assay` when the state records the
  ITSD outcome. Both kept, KIT-242's inside the list, KIT-255's after it.
- NEWS auto-merged, staging's bullets first. KIT-255's interim sentence ("Where a sample's values were missing altogether
  in one assay, that assay may not have been normalised between samples (see KIT-242; until KIT-242's first part is
  merged).") is removed: with Part A such an assay is normalised on its other samples, so the sentence is now untrue.
- Ledger: kit-242 adds no ledger row, so staging's rows stand as they are.

**R2's base.** Branch `kit-242-r2base` = 9bb5ef0d: afb27738 with Part B's code commits (5c856dff, the R part of
4c8e3d87) reverted; tests identical to the head's (`git diff` lists only `R/` and `NEWS.md`). Head: afb27738.

**Files: 190** (`r2-files.txt`, recounted on the merged tree): the 186 of P0, plus the 6 that the scope's grep now finds
(test-kit213-study-parameters-labels.R, test-kit221-da-several-comparisons.R, test-kit255-omics-correlation-assays.R,
test-kit259b-second-correction.R, test-kit264-trx-export-current.R, test-kit271-limpa-table-name.R), minus the two
runner-library files. 12 shards per side, check on both sides; `env` `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` and
`MULTISCHOLAR_KIT105_W3_DUMP=1`, `extra_packages` limpa, mixOmics, sva, `file_timeouts` 120 for test-prot-norm-core.R
and the GUI-entry files (as train 19). Line: 59.33 expected / 120.0 ceiling core-hours.

**Predictions, restated on the merged tree** (the table and the `test-omic-norm-module.R` keys above stand, with):
- `test-omic-norm-module.R`: on staging 631/8 (KIT-255's eight). Base fails the union of those 8 and Part A's 56 (the
  count is observed on base: 56 if the 8 are among Part A's keys, up to 64 otherwise); head fails the same set. Verdict
  message-changed, on exactly the log, provenance and outputs keys listed above (and the two `reset outputs` keys not
  predicted); every other failing key, KIT-255's included, has the identical message on both sides.
- `test-prot-norm-core.R`: 1547/307 on staging; base and head identical (Part C is on both sides). Same.
- `test-kit236-report-sentences.R` (293/0) and `test-metab-s4-qc-shared.R` (35/0): pass on both sides. Same.
- New files: `test-kit255-omics-correlation-assays.R` (360/0: its legacy RUV is Manual on complete objects, every assay
  corrected, so the head's result equals the base's), `test-kit271-limpa-table-name.R` (78/0), `test-kit221-…` (124/0),
  `test-kit259b-…` (90/0; its `ruv_applied` is a stub report value), `test-kit213-…` (56/0; transcriptomics),
  `test-kit264-trx-export-current.R` (transcriptomics; no lipid or metab RUV): same on both sides.
- Any other verdict stops the run.

### R2 (run 45738934354, base 9bb5ef0d, head afb27738; 1102.3 job-minutes, 73.49 core-hours; ceiling 120.0): observed

- **File verdicts as predicted.** fixed: `test-kit242-normalisation-fallback.R` (base 370/84/4 errors, head 526/0),
  `test-workflow-omics-worker.R` (base 4 failures, head 0), `test-module-ci-metab-norm.R` (base 4, head 0).
  message-changed: `test-omic-norm-module.R`, `test-workflow-omics-declarations.R` (KIT-073),
  `test-workflow-omics-binding-steps.R` (KIT-094). The other 184 files same; check identical (6 findings on both).
  `test-omic-norm-module.R` 583/56 on both sides (KIT-255's eight are among Part A's 56 keys).
- **Expectation level, `test-omic-norm-module.R`: as predicted** (4 log, 12 provenance, 10 outputs keys changed). The two
  keys not predicted, manual `reset outputs` (lipid, metab), changed: the log survives the reset.
- **KIT-094: as predicted**, the `records` expectation changed; `report_sources` identical.
- **KIT-073: one departure from the written expectation list.** The two `values` expectations changed as predicted; the
  `registry` expectation's message is identical on both sides, where a change was predicted. Inferred cause: the
  registry's difference display is cut after its first differences (it begins with
  `lipid_between_sample_normalisation`, `lipid_duplicate_resolution`, `lipid_intensity_filter`), so the RUV
  declarations' new `fixed` row is never printed; the actual value is expected to differ but the message cannot show it.
  Not verified. Stopped and reported under the stop rule.
- The E25m row fails by value on the base, as the change intended; the two declared edits fail on the base and pass on
  the head.
- Used: 148.31 + 73.49 = 221.80 of 901.33 core-hours.

## After R2: the coordinator's rulings (2026-11-26), no run

**(1) KIT-073's `registry` message: a prediction error of the message, verdict unaffected.** Established by reading
(the failure texts of run 45738934354 and the code):
- The expectation is `expect_identical(current$registry, frozen$registry)` (`tests/testthat/test-workflow-omics-declarations.R:45`),
  against the stored dump `tests/testthat/fixtures/kit073-omics-declarations-frozen.rds` (`:7-10`). The failure text is
  printed in full (6,664 characters, identical on base and head in `failures.txt` and `failures-base.txt`); it is not
  cut. It ends with "`actual$metabolomics$metab_ruv$declaration$fixed` is a list / `expected$...$fixed` is NULL", and
  reports `lipid_ruv`'s the same way: the stored dump predates every `fixed` list on the RUV declarations (KIT-224's
  first one), so the comparison stops at "list against NULL" for that entry and never prints the list's rows.
- The value does differ: head `fixed = list(ruv_feature_levels = "restored", ruv_control_grouping = "kept", ruv_k_bound =
  "sample_size", ruv_fallback = "by_cause")` (`R/utils_workflow_omics_declarations.R:790-791`; the declared design's
  `fixed` gains `ruv_fallback` at `:1102`); base 9bb5ef0d the same list without `ruv_fallback`
  (`R/utils_workflow_omics_declarations.R:760` there). So `current$registry` differs between base and head; the message
  cannot show it because an earlier difference of the same entry is what is printed. My earlier inferred cause (a display
  cut after the first differences) was wrong; this is the same situation as KIT-258's (VAL-338).

**(2) The two keys declared "not predicted" (manual `reset outputs`, lipid and metab) changed.** What it means: the
tab's log (the `norm_log` output) still holds the RUV run's closing line after the reset, so the reset does not clear
the log. The ticket's text says nothing about the reset or the log; listed for the auditor (whether a reset should clear
the log is outside KIT-242; its behaviour here is the existing one, only the line's text changed).

**(3) Run order (the coordinator's ruling).** R3 is not sent: the implementation audit comes first, and A5's full suite is
met by the merge train's run (KIT-242 with KIT-260, the full suite on the merged tree), so that one full suite is not
run twice.

**NEWS clause for the train (prepared, not committed here; KIT-268's bullet is on staging since train 32, not on this
branch).** Staging's text "(for lipidomics and metabolomics once their analysis record names the method per comparison,
KIT-242) (KIT-268)" becomes "; the lipidomics and metabolomics analysis record now names the method per comparison and
assay (KIT-268, KIT-242)". KIT-236's bullet ("Not yet changed (they wait for KIT-242): the sentence naming the
between-sample normalisation method applied, the internal-standard sentence and the method named in the Figure 1
caption") stays true: KIT-242 writes the record those sentences will read but changes no report template (section 3,
out of scope, KIT-236).

**Ledger.** A6 owes a row: `VAL-<merge>` written after VAL-358 in `docs/methods-decision-ledger.md` (the decisions as
implemented, the evidence runs, the precedents, the classification and what is not measured).

### A3 (Part B's methods)

`git grep -n "log_warn(" -- R/func_lipid_s4_ruv_methods.R R/func_metab_s4_ruv_methods.R` at 8162ea52 (55 each, twins):
- `ruvIII_C_Varying`: under the current rule the method hands every call to `omicsRuvCorrectByCause`, which runs the
  method body under the earlier rule only for the assays it corrects; every skip there (`return(NULL)` after a
  `log_warn`) leaves the assay absent or unchanged, which `omicsRuvCorrectByCause` turns into `omics_ruv_failed`
  (`correction_error`), and the invalid-k, malformed-assay and no-control cases are decided before the body runs. Two
  logged-only removals remain, as decided: features with values returned all non-finite (`:713-715` area, G9) and
  samples returned all non-finite (`:719-722` area, adjacent item 8, not decided; an excluded empty sample never reaches
  it, it is put back as missing).
- `ruvCancor`: a skipped assay returns no plot; under the current rule the search records it as `error_cancor`, a failed
  computation when no percentage is valid (stop). The control alignment at `:219` no longer skips (Q6.1).

## Acceptance (ticket section 9)

| # | Criterion | Evidence | State |
|---|---|---|---|
| A1 | G0 observes C1 to C14 as predicted (or the ticket is revised) | `g0-*` on `kit-242-g0`; r1.3 | met (C4 revised the ticket, r1.3) |
| A2 | Every E row fails at base, G rows pass, all pass on the fix | R1a, R1b-A, R1c, R1b-B, R1d; R2 base/head | met (R2: the new file fixed, 526/0 on the head) |
| A3 | No warning followed by an unhandled skip in the four S4 methods | the A3 sections above | met on the application's path, with the two logged-only removals named; except the script path of audit I07: `omicsRuvCorrectByCause` hands a k or control list that does not name every assay to the earlier rule, which drops that assay silently |
| A4 | P0 lists the moving expectations; R2 as 8.2 | P0 45723354055 (`moved.tsv`); R2 45738934354 | not met for the 18 expectations of audit I01, which move outside the owner-approved files (awaiting the owner, asked 2026-11-26); otherwise met, with one message-level prediction error (KIT-073 `registry`), verdict unaffected |
| A5 | R3 full suite: no difference other than A4's | by the merge train's run, coordinator's ruling | owed to the train |
| A6 | NEWS, ledger row, `outcome.md` with "For the change report"; `man/`, `NAMESPACE`, `DESCRIPTION` unchanged | NEWS (two KIT-242 bullets), ledger `VAL-<merge>`, this file | met except that OD-4 (c)'s note wherever an older record is read is not implemented (audit I03), and `VAL-<merge>` now says so; `git diff aa119308 -- man NAMESPACE DESCRIPTION` empty |
| A7 | Hand-offs recorded by the coordinator | `DISPATCH-STATE.md` | the coordinator's |

## Stored values moved, by file (for the stored-reference rows; VAL-303, VAL-352)

No stored value is edited beyond the two VAL-303 pins and the three VAL-352 edits. Failing until the end-of-project
re-record (VAL-303, VAL-352 (2) to (4)):
- `tests/testthat/test-omic-norm-module.R` (frozen captures `:553`-`:1539`, compared at `:1549`, call sites
  `:1560`-`:1585`): 56 expectations: 22 provenance, 12 states and 6 final-PCA digests (whole S4 states, moved by the
  per-assay ITSD arguments), 12 outputs and 4 log lines. 12 are VAL-303 (1)'s, the other 44 VAL-352 (4)'s.
  Part B also moves, inside those 56: log (manual run, automatic run), provenance (manual run, skip_correlation, reset;
  automatic run, final_plot, final_plot_failing), outputs (manual run, skip_correlation, reset; automatic run,
  final_plot, final_plot_failing), lipid and metab (observed R2).
- `tests/testthat/test-workflow-omics-declarations.R:45` and the `values` comparisons (KIT-073, fixture
  `fixtures/kit073-omics-declarations-frozen.rds`): registry (four marker rows) and the RUV values (five new fields);
  VAL-352 (3).
- `tests/testthat/test-workflow-omics-binding-steps.R:197` (KIT-094, fixture `fixtures/kit094-characterisation.rds`):
  the `records` part (the marker rows; the RUV step's per-assay outcomes); VAL-352 (3).
- GUI-entry lanes (`tests/testdata/gui-entry/lipid_canonical.json`, `metab_lc.json`, `metab_gc.json`,
  `metab_combined.json`; `test-gui-entry-lipidomics.R`, `test-gui-entry-metabolomics.R`): the internal-standard log and
  record captures (Part A); no RUV in these lanes; VAL-352 (2).
Awaiting the owner (asked 2026-11-26; audit I01): 18 expectations in three files that pass on staging and that Part A
turns red; no owner row covers them (not VAL-303, not VAL-352, not the standing re-record rule, which covers only
references already failing):
- `tests/testthat/test-workflow-lipid-declarations.R:138`, `:188` and `test-workflow-metab-declarations.R:138`, `:188`
  (2 each): `itsd_applied` TRUE on a fixture whose ITSD matches no standard, and in the report overlay.
- `tests/testthat/test-kit172-downstream.R` T10 (`:189-238`, 14): the tab's warning and notification counts and the
  first warning's text (the fixture matches no standard, so the tab also warns that ITSD was not applied).
Edited, with authority: `test-workflow-lipid-declarations.R:47`, `test-workflow-metab-declarations.R:47` (VAL-303 (2));
`test-metab-02w-…:226-240` (VAL-352 (1)); `test-workflow-omics-worker.R` RUV branch and `test-module-ci-metab-norm.R:200-204`
(VAL-352 (5)); `test-workflow-omics-worker.R:15-23` (Part A's declared edit, r1.6).

## Not checked (final)

- The full suite on the merged tree (A5): by the merge train's run.
- OD-4 (c)'s note wherever an older record is read (record viewer, report overlay, KIT-236's sentences): not implemented
  here; the markers are (audit I03; the owner told).
- Audit I04: on the proteomics `normalization-ruv` route an empty sample leaves the design with no notice and no record
  (the warning `omics_normalisation_sample_removed` is not captured by `diaNormalizationOperation`); owed by r1.6, not
  done; awaiting the owner's decision (asked 2026-11-26).
- Audit I05: the study-parameters file's ITSD lines are not per assay (owed by r1.6 "with Part B"; not done).
- Audit I06: the declared design records the excluded sample only in its factor records (no record outcome, no notice;
  its warnings escape uncaptured); the gate is off.
- Audit I11: no test for `protein_normalisation_fallback` in the two proteomics declarations, its earlier-token scale
  case or a pre-KIT-242 proteomics record's replay; E18 to E23 under the earlier RUV token (only E14 is shown); no row for
  the notice naming a one-value sample; `test-omics-module-globals.R` not in R2's scope.
- The declared design through the tab and the managed worker with an empty sample (E20 is kernel-level).
- RUV-III-C's own behaviour on a sample left out and put back (item iii's "not checked"): the put-back column is
  all-missing by construction; the DA on it is limma's per-feature fit.
- Whether any expectation already failing on staging hides a numeric change: R2's per-expectation comparison shows only
  the listed message changes.
- Real data (KIT-188 and the KC sets are not used); the frequency of empty samples and of no-elbow searches in real data.
- Whether a reset should clear the tab's log (listed for the auditor above).

## For the change report (final)

- Lipidomics and metabolomics between-sample normalisation left an assay unnormalised when one sample had no value or the
  method failed, while the tab said complete and the session summary said it was done (legacy, on `main`); the analysis
  record saying "applied" is new in 0.9 (legacy code, one trigger new in 0.9: IMP-017). Now the sample is left out of that assay only (removed from the design when empty everywhere) and failures
  stop. Evidence: G0 C1-C3, R2a-2, R2.
- ITSD with no matched standard (a likely default case) normalised nothing for that assay and said "applied"; now
  recorded per assay with the standards used, shown in the tab and the session summary. Evidence: G0 C3, C10.
- Legacy RUV-III-C dropped an assay it could not correct, or returned the input when every assay failed (legacy, on
  `main`); the analysis record's `ruv_applied` TRUE for it is new in 0.9. The control alignment failed on any feature
  the control selection did not list (metab since 2025-05-06, lipid since 2026-01-18). Now by cause (VAL-242 Q6.1, Q6.2 (c); VAL-290). Evidence: G0 C4, C12,
  R1b-B, R1d, R2.
- The declared RUV design refused an assay with an empty sample; now computes its factors on the other samples and
  leaves the sample out of that assay's DA. Evidence: E20 (R1d, R2).
- Proteomics Scale normalisation with an empty sample made the table all-missing; the sample now leaves the design
  (Part C). Evidence: G0 C9, E17.
- Every change has a replay marker, so records made before replay as made.

## After the implementation audit (AUDIT-KIT-242-impl-20261126, I01 to I11): record fixes; no code, no test change

- R2's head (afb27738) is the tip's code: every later commit touches only `docs/`. The merge afb27738 resolved **two** R
  files by hand (I08), not one as said in an earlier hand-back: `R/utils_workflow_analysis_record.R` (the marker list)
  and `R/utils_workflow_saved_exports.R` (the session list).
- P0's class "already failing on staging" is 951 (926 + 25) and P0 cost 35.95 core-hours (539.2 x 4 / 60), corrected above
  (I08). Ticket total unchanged to rounding (221.81).
- Corrected above: the 18 expectations of I01 (awaiting the owner, asked 2026-11-26; not on any re-record list),
  acceptance A3 (I07), A4 (I01) and A6 (I03), "Not checked" (I03's note, I04, I05, I06, I11), "For the change report"
  (I08), the description of the 56 digests (I09 (a)).
- The lanes: no RUV marker enters them (they skip RUV, `tests/testthat/helper-gui-entry.R:615`); VAL-352 (2)'s mention of
  "the RUV marker" should not be repeated in the stored-reference row (I09 (c)).

### Prepared for the owner (text only; not committed as code or tests)

**(i) I01 under option (a): the 18 edits, old and new.**

`tests/testthat/test-kit172-downstream.R`, first T10 test (`:207-209`, per route; 8 expectations over lipid and metab, in
process and managed). Old:

    warnings <- Filter(function(note) identical(note$type, "warning"), run$notes)
    expect_length(warnings, 1L)
    expect_identical(warnings[[1L]]$message, text, info = route)

New:

    # KIT-242 (VAL-242 OD-2 (b); owner <date>, VAL-<n>): the fixture matches no internal standard, so the tab also warns
    # that ITSD was not applied; T10 counts the scale warning by its own text.
    warnings <- Filter(function(note) identical(note$type, "warning") && identical(note$message, text), run$notes)
    expect_length(warnings, 1L)
    expect_identical(warnings[[1L]]$message, text, info = route)

Second T10 test (`:236`, three cases per class; 6 expectations). Old:

    expect_length(Filter(function(note) identical(note$type, "warning"), run$notes), 0L)

New:

    # KIT-242 (owner <date>, VAL-<n>): no scale warning; the ITSD notice for this fixture is not counted here.
    expect_length(Filter(function(note) identical(note$type, "warning") &&
      grepl("Scale normalisation rescales", note$message, fixed = TRUE), run$notes), 0L)

`tests/testthat/test-workflow-lipid-declarations.R:138` and `test-workflow-metab-declarations.R:138`. Old
`expect_true(by_state$lipid_itsd_norm$outcomes$itsd_applied)` (metab: `metab_itsd_norm`); new
`expect_false(by_state$lipid_itsd_norm$outcomes$itsd_applied)  # KIT-242 (owner <date>, VAL-<n>): no standard matched`.
`:188` in both: old `expect_true(overlay$itsd_applied)`; new `expect_false(overlay$itsd_applied)  # KIT-242: no standard
matched, so the overlay says not applied`. Shown failing on staging and passing on the train's head.

**(ii) I04, the smallest repair.** By reading, the in-process twin `runProtNormBetweenSamplesStep`
(`R/mod_prot_norm_workflow_helpers.R:15-30`) is no longer reached by the application: `grep` finds no caller in `R/`,
only two characterisation tests (`test-prot-05b-norm-module-contracts.R:1006-1029`, `test-prot-05d-norm-module-characterization.R:1252-1262`).
So only `diaNormalizationOperation` needs it:
- `R/utils_workflow_dia_products.R:387`: wrap the call as the limpa route does (`R/utils_workflow_prot_declarations.R:1252-1257`):

      removed <- character()
      normalized <- withCallingHandlers(normaliseBetweenSamplesFn(object, normalisation_method = parameters$norm_method),
          omics_normalisation_sample_removed = function(w) {
              removed <<- unique(c(removed, w$sample))
              invokeRestart("muffleWarning")
          })

- `:510`: `if (length(removed)) value$samples_removed <- removed` (added only when non-empty, so other products are
  unchanged).
- The record (`R/utils_workflow_analysis_record.R:155-170`, the `normalization-ruv` outcomes): `if
  (length(saved$samples_removed)) outcomes$samples_removed <- saved$samples_removed`; the tab's step text and a warning
  notification "Warning: samples <samples> have no value and were removed from the design." (the limpa route's text,
  `diaProteinNormalisationText`), in both runners' result handlers for this step.
- Test (proteomics, E17 through the declaration): a protein object with B_6 empty, `runSavedOperation` of
  `diaNormalizationOperation` with `norm_method = "cyclicloess"`, `ruv_mode = "skip"`: the value's `samples_removed` is
  "B_6", the record's outcomes name it, the design lacks it. Fails on the base by value (`samples_removed` NULL);
  `test-prot-norm-core.R` stays `same` (no fixture holds an empty sample).

**(iii) I02, the two candidate behaviours (lipid and metab; under the current rule only).** Common part:
- The control selection records why it skipped an assay, under the current rule: lipid
  `R/func_lipid_s4_ruv_control_methods.R:187-191` (groups with fewer than two samples), `:148-162` (duplicated ids),
  `:167-172` (too little data); metab `R/func_metab_s4_design_control_methods.R:270-273` and twins. Sketch: return the
  skipped assays' reasons as an attribute of the list (`attr(final_control_list, "skipped") <- list(A2 = list(reason =
  "group_counts", groups = c(B = 1L)))`), set at `:348` (metab `:431`).
- Automatic search: lipid `R/func_lipid_norm_ruv_optimization_helpers.R:191-203`, metab `R/func_metab_norm.R:714-726`:
  when the control list is named and lacks the assay, return a row with status `skipped_no_controls_selected` and the
  attribute's reason (instead of taking the whole list, which makes `sum()` fail at `:203` / `:726`).
- Manual search: lipid `:421-425`, metab `:944-948`: the same case returns `success = FALSE` with `skip_reason`.
- `omicsRuvSearchOutcomes` (`R/utils_workflow_omics_declarations.R:804-845`): map `skipped_no_controls_selected` (and a
  Manual result with `skip_reason`) by the reason; duplicated ids keep today's stop with "Resolve Duplicates".

(a) Stop with the true cause:

    if (identical(reason$reason, "group_counts")) return(list(action = "stop", reason = "group_too_small",
        text = sprintf(paste("RUV-III was not run: in assay %s the groups %s have fewer than two samples with values, so",
            "control features cannot be chosen. Remove the sample, choose another grouping, or run without RUV."),
            assay_name, paste(names(reason$groups), collapse = ", "))))

(b) Keep that assay without RUV, recorded as `no_usable_controls`:

    if (identical(reason$reason, "group_counts")) return(list(action = "keep", reason = "no_usable_controls",
        text = sprintf(paste("RUV-III was not applied to assay %s: the groups %s have fewer than two samples with values",
            "there, so control features cannot be chosen. The assay keeps its values without RUV."), assay_name,
            paste(names(reason$groups), collapse = ", "))))

Test rows (lipid and metab, Automatic and Manual; E27): A2's group B left with one sample holding values
(`v$A2[, c("B_2", "B_3", "B_4", "B_5", "B_6")] <- NA`), the chain to the normalised state, then legacy RUV. Base (the
tip): Automatic stops with R's "invalid 'type' (list) of argument", Manual with `search_error` and the wrong cause. Head
under (a): `omics_ruv_failed`, `assay` "A2", `reason` "group_too_small", the message naming group B; under (b): the step
completes, A2 equal to its input, A1 corrected, `ruv_not_applied_reason_per_assay$A2` "no_usable_controls", the notice
naming group B. Plus KIT-275's VE shape (group B with one sample in every assay): under (a) a stop naming every assay;
under (b) every assay kept without RUV, `ruv_applied` FALSE.

## VAL-373: the owner's answers and the one bounded correction (Q-C)

VAL-373 (origin/v0.9-staging, 2026-11-27): Q-A (b), Q-B (a), Q-C one bounded correction for two things only (the DIA
normalisation route's notice and record of a removed sample; Q-A's behaviour on both search paths), Q-D the OD-4 (c) note
is a follow-up.

**Tests-only commit a8c299c2** (nothing under `R/`): E17 through `diaNormalizationOperation` (kernel, commit outcomes,
the tab's notice), E27 (lipid and metab, Automatic and Manual; seed 242027) and KIT-275's VE shape (lipid and metab;
seed 242227), and the 18 edits of Q-B exactly as prepared (comments cite VAL-373).

### Predictions for R1e, written before the run (a8c299c2: the code of 11d3150a, unchanged since afb27738, plus these tests)

Files: `test-kit242-normalisation-fallback.R`, `test-kit172-downstream.R`, `test-workflow-lipid-declarations.R`,
`test-workflow-metab-declarations.R`; 1 job, ceiling 60 job-minutes.
- `test-kit242-normalisation-fallback.R`: every earlier row passes (526 as on R2's head). E17 fails by value on four
  expectations (the removed-sample warning escapes the kernel: `kept` length 1; `samples_removed` NULL; the commit's
  `outcomes$samples_removed` NULL; no warning notification naming B_6). E27 Automatic (lipid, metab) fails at
  `expect_false(inherits(run, "error"))`: the step stops with R's "invalid 'type' (list) of argument". E27 Manual (lipid,
  metab) fails there too: the step stops with `omics_ruv_failed` (`search_error`, "could not read the assay"). The VE
  shape (lipid, metab) fails there with the Automatic R error. No error is raised by any row (all failures are by value).
- `test-kit172-downstream.R`: passes in full (84 passed, 0 failed): T10 now counts the scale warning by its own text,
  which the code already gives once; the ITSD notice is not counted.
- `test-workflow-lipid-declarations.R`, `test-workflow-metab-declarations.R`: pass in full (259 passed, 0 failed each): on
  this code `itsd_applied` is FALSE on the fixture whose ITSD matches nothing (`:138`) and in the overlay (`:188`); `:47`
  as re-pinned.
- A miss is a result: recorded and handed back.
