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
