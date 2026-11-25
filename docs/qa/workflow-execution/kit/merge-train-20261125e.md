# Merge train 5 of 2026-11-25: KIT-239, KIT-214, KIT-226, KIT-241

Coordinator's instruction after train 4 landed (staging e0747fc1, ledger rows through VAL-233): same
procedure as train 3 (`merge-train-20261125c.md`).

## Train

Base: `e0747fc198817596e3e66397a87e06b2f869c6a7` (`origin/v0.9-staging`).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-239 | `origin/kit-239` f41ee16b | 72c7e708 | `VAL-255` | AUDIT-KIT-239-impl-20261125, closure check |
| 2 | KIT-214 | `origin/kit-214` 13677193 | 032abdd7 | `VAL-256` | AUDIT-KIT-214-impl-20261125 |
| 3 | KIT-226 | `origin/kit-226` 05c4376a | b4422d1c | `VAL-257` | AUDIT-KIT-226-impl-20261125, closure check |
| 4 | KIT-241 | `origin/kit-241` 2eb80db9 | cc3d4bdc | `VAL-258` | AUDIT-KIT-241-impl-20261125 |

Head (bound by both runs): `cc3d4bdc` (full hash in the requests). This file and the records commit
after the runs are under `docs/qa/` only.

Conflicts met:
- `docs/methods-decision-ledger.md` in all four merges: staging's rows through VAL-233 kept; each
  member's row after them as `VAL-<train5-k>` with a clean status cell. KIT-239's placeholder also
  replaced in INT-28's implementation column (closure check, table 3). Each member's `outcome.md`
  placeholder replaced likewise.
- `NEWS.md` in merge 4 only (KIT-241's Fixes entry against the entries already on staging): every
  entry kept, under the one unreleased heading. No NEWS entry mentions KIT-222's merged-samples
  defect.
- No conflict under `R/`, `tests/`, `inst/`, `man/`, `DESCRIPTION` or `NAMESPACE`;
  `R/mod_omic_enrich_panels.R` merged automatically between KIT-226 and KIT-241. `DESCRIPTION`
  gains `zip` in Imports (KIT-241); its Version line is unchanged.

Path lengths (checked by script before sending): 82 paths are new since the base; the longest is 89
bytes with the `MultiScholaR/` prefix; none exceeds 100.

## Declaration 1: run-tests (written before it is sent)

- Purpose: each member's check on the merged tree, done once on the train.
- Scope: the union of the members' merge-check lists, 92 files, base and head, R CMD check on both
  sides, `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` on both sides:
  - KIT-239: the 19 files of its R2 (`kit239-r2-09286c2d.json`);
  - KIT-214: the 49 files of its R2 (`kit-214/r2/compare.tsv`);
  - KIT-226: the 23 files of its R2 (`kit226-r2-77ab881f.json`), which include the new file,
    `test-omic-enrich-module.R`, `test-omics-enrich-core.R` and `test-gui-entry-proteomics.R`;
  - KIT-241: its register's list (section 6): `test-kit241-enrich-download-zip.R`,
    `test-kit215-not-tested-lists.R`, `test-enrichment-significance.R`, `test-omic-enrich-module.R`,
    `test-module-ci-prot-enrich.R`, `test-module-ci-manifest.R`, the two `test-prot-11` files,
    `test-e2e-proteomics-enrichment.R`, `test-kit226-gprofiler-columns.R`.
  8 shards per side; per-file limits of 120 minutes for the three GUI-entry files and
  `test-workflow-partitions.R` (kill limits, as the members' R2s).
- Ceiling: 2,020 job-minutes, the sum of the members' merge-check ceilings: KIT-239 420 (its R2,
  the check its closure names); KIT-214 600 (R2b, ticket section 10); KIT-226 540 (R2, whose scope
  the closure check names); KIT-241 460 (R2).
- Predictions (base values from train 4's compare on the landed tree):
  - KIT-239: `test-integration-results.R` 169 to 258 passed; `test-integration-record.R` 133 to
    167; `test-integration-gui.R` 397 to 401 (all `count-changed`, 0 failed);
    `test-integration-model.R` 263 `same`; `test-integration-removed-factors.R` 520 `same`, 3
    skipped; the other 14 files `same` (the MOFA2 files skip on both sides).
  - KIT-214: `test-kit214-sample-columns.R` only on head, 176 passed, 0 failed;
    `test-workflow-omics-declarations.R` 45/23 and `test-workflow-omics-binding-steps.R` 58/1,
    counts unchanged, `message-changed` by exactly KIT-214's lines (declarations: `declaration$fixed`
    a list instead of NULL for `lipid_duplicate_resolution` and `metab_duplicate_resolution`;
    binding steps: `steps[[4]]$parameters` length 1 instead of 0 in both omics; no other train member
    adds a marker); `test-workflow-kit-harness.R` 734/1 and `test-workflow-operation-registry.R`
    1095/2 `same`; every other file `same`.
  - KIT-226: `test-kit226-gprofiler-columns.R` only on head, 38 passed, 0 failed;
    `test-omic-enrich-module.R` exactly the base's failures (1051/19: clusterprofiler 8,
    kernel_clusterprofiler 9, KIT-113 W1 2; none in gprofiler, kernel_gprofiler, organism_filter,
    restore, A107-04 or A107-08), `same`; `test-omics-enrich-core.R` the base's list only, `same`;
    `test-gui-entry-proteomics.R` 189/4 `same` (`prot_dia` back to its four integration-file
    failures); `test-enrichment-significance.R` passing (one assertion replaced, count unchanged);
    the other R2 files `same`. The g:Profiler W1 download leaves at the re-frozen values
    (summary `c5b5570d…`, `gprofiler2_results.tsv` `71e8ee7c…`) pass under `zip::zipr`.
  - KIT-241: `test-kit241-enrich-download-zip.R` only on head, 15 passed, 0 failed; every other
    file of its list `same`; `test-e2e-proteomics-enrichment.R` skipped on both sides.
  - R CMD check: findings identical to the base's (or fewer on head); `zip` never in an
    "'::' or ':::' import not declared" line; the portable-file-names NOTE lists the same number of
    paths on both sides.
- Decision: every prediction holds: that member's part is green. A member whose prediction fails is
  named with files and counts, not repaired here.

## Declaration 2: the MOFA lane for KIT-239 (written before it is sent)

- Purpose: KIT-239 closure check, "What the train must hold" item 2 (train 4 did not run it).
- Request: `mofa-package.yml` on branch `run-mofa-smoke-train5` (a copy of `run-mofa-smoke-k239`,
  the way KIT-239's R3 45577930225 was requested: the workflow's trigger branch set to the copy,
  then `.github/mofa-smoke/package-params.json`): `package_commit` the train's head,
  `tests` `test-integration-workbook-mofa.R` and `test-integration-record-mofa.R`, `record` false,
  `runs` 1.
- Ceiling: 60 job-minutes (KIT-239's R3 row).
- Predictions:
  - `test-integration-workbook-mofa.R`: 3 tests, 49 expectations passed, 0 failed, 0 skipped. Read
    at the base and at the head: KIT-234 archived the workbooks, but this file stayed. It holds three
    tests: MOFA2's defaults; the module reproducing the reference workbook's configuration on
    synthetic data; and the ported "KIT-143 (KIT-234): the twelve integration functions run end to
    end on the corpus with MOFA2 and repeat the model from the record", which renders the report.
    The file is unchanged since KIT-234's L run 45524568014 at 652c962b (git log), which gave 49
    passed and 0 skipped. KIT-234's register I07 names this run as the train's to hold, because
    KIT-239 changes the report template and `R/utils_integration_record.R` it renders.
  - `test-integration-record-mofa.R`: 29 passed, 0 failed, 0 skipped, with the line
    "KIT-239 L8: likelihoods gaussian, gaussian, gaussian, gaussian, gaussian, gaussian; groups 1".
  - Any skip is a stop (KIT-239 stop rule 7).

## Result

Both runs were sent once.

### run-tests: 45598541744

The request was `merge-train-20261125e-cc3d4bdc`, run-tests commit 9f083f47, head cc3d4bdc against base e0747fc1.
- Conclusion: success.
- 578.1 job-minutes (check 30.1).
- Tables in `merge-train-20261125e/`: `compare.tsv`, `tests.tsv`, `tests-base.tsv`, `check-compare.txt`,
  `check-summary.txt`, `meta.txt`.

Observed: of 92 files, 84 are `same` and 8 are not.

| file | base p/f/e/s | head p/f/e/s | verdict | member, prediction |
|---|---|---|---|---|
| test-integration-results.R | 169/0/0/0 | 258/0/0/0 | count-changed | KIT-239, 169 to 258: held |
| test-integration-record.R | 133/0/0/0 | 167/0/0/0 | count-changed | KIT-239, 133 to 167: held |
| test-integration-gui.R | 397/0/0/0 | 401/0/0/0 | count-changed | KIT-239, 397 to 401: held |
| test-kit214-sample-columns.R | absent | 176/0/0/0 | only-on-head | KIT-214, 176/0: held |
| test-workflow-omics-declarations.R | 45/23/0/0 | 45/23/0/0 | message-changed | KIT-214: held (below) |
| test-workflow-omics-binding-steps.R | 58/1/0/0 | 58/1/0/0 | message-changed | KIT-214: held (below) |
| test-kit226-gprofiler-columns.R | absent | 38/0/0/0 | only-on-head | KIT-226, 38/0: held |
| test-kit241-enrich-download-zip.R | absent | 15/0/0/0 | only-on-head | KIT-241, 15 passed: held |

Predicted `same` and observed `same` (passed/failed/skipped, on both sides):

| File | Result |
|---|---|
| `test-integration-model.R` | 263/0 |
| `test-integration-removed-factors.R` | 520/0/3 |
| `test-integration-workbook-mofa.R` | 0/0/3, skipped in run-tests |
| `test-integration-record-mofa.R` | 0/0/1, skipped in run-tests |
| `test-workflow-kit-harness.R` | 734/1 |
| `test-workflow-operation-registry.R` | 1095/2 |
| `test-omic-enrich-module.R` | 1051/19 |
| `test-omics-enrich-core.R` | 232/9/1 |
| `test-gui-entry-proteomics.R` | 189/4 |
| `test-enrichment-significance.R` | 181/0 |
| `test-e2e-proteomics-enrichment.R` | 0/0/2 |
| `test-kit215-not-tested-lists.R` | `same` |

In `test-omic-enrich-module.R` the head's 19 failures are clusterprofiler 8, kernel_clusterprofiler 9 and KIT-113 W1 2, the
same as the base's. None is in gprofiler, so the re-frozen g:Profiler download leaves (`c5b5570d…`, `71e8ee7c…`) pass
under `zip::zipr`.

KIT-214's message changes, failure records base against head (temporary paths masked):
- In the declarations file, the head adds exactly 6 lines: `declaration$fixed` is a list instead of NULL for
  `lipid_duplicate_resolution` and `metab_duplicate_resolution`.
- In the binding-steps file, the head adds exactly 12 lines: `steps[[4]]$parameters` is length 1 instead of 0, and its
  `[[1]]` is a list instead of absent, in both omics. These lines are identical to KIT-214's own R2 record
  (`kit-214/r2/compare-failures-binding-steps.txt`, lines 51-55 and 71-75).
- Nothing else moved.

R CMD check:
- 0 ERRORs, 3 WARNINGs, 3 NOTEs on both sides; all 6 findings identical, none only on head.
- The "'::' or ':::' import not declared" line names only `cli`, on both sides; `zip` is never named.
- The portable-file-names NOTE lists 238 paths on each side.

### MOFA lane: 45598538505

The branch was `run-mofa-smoke-train5` (56154c47), package cc3d4bdc.
- Conclusion: success; exit 0 after 307 s.
- 14.6 job-minutes of the 60 ceiling.
- Tables in `merge-train-20261125e/mofa/`.

Results:
- `test-integration-workbook-mofa.R`: 3 tests, 8 + 29 + 12 = 49 passed, 0 failed, 0 skipped. Held.
- `test-integration-record-mofa.R`: 29 passed, 0 failed, 0 skipped. The console line reads "KIT-239 L8: likelihoods
  gaussian, gaussian, gaussian, gaussian, gaussian, gaussian; groups 1". Held.

### Verdict

Every prediction of KIT-239, KIT-214, KIT-226 and KIT-241 held, and no other file differs. The check is identical. TRAIN
GREEN.

Job-minutes: 578.1 against the 2,020 ceiling, plus 14.6 against the MOFA lane's 60.

### Staging and numbering

While the runs were out, staging moved by ledger-only commits: e0747fc1 → 29065062, which touches only
`docs/methods-decision-ledger.md`, with more such commits announced. The train's base therefore differs from staging in
the ledger alone, and the check stands. The train lands by a merge. Its four rows stay `VAL-255` to
`VAL-258`, numbered by the coordinator at landing, in train order.

### Not checked

- The trees after the runs' head: docs only.
- Files outside the 92.
- The members' shared full suites, which are owed after the train (KIT-226 section 6; KIT-241 A4).
- The merge with the current staging.

## Merge record (coordinator)

Landed on `v0.9-staging` on 2026-11-25 by a merge of this branch at 5cf10ded into staging 7327607f.
Staging had moved from the train's base e0747fc1 by two commits that touch only the methods ledger
(the owner's decisions of that evening, rows VAL-234 to VAL-254). The only conflict was the ledger:
staging's rows were kept through VAL-254 and the train's rows numbered after them in train order:
VAL-255 (KIT-239), VAL-256 (KIT-214), VAL-257 (KIT-226), VAL-258 (KIT-241); the numbers were written
into the ledger (KIT-239's also in INT-28's implementation column), each member's `outcome.md` and
this file. One coordinator's row was added at landing, VAL-259: KIT-214's marker shows in the two
laptop-recorded snapshots of VAL-181 and VAL-189, which are not edited and wait for the single
end-of-project re-record. The landed tree differs from the train's tip only in the ledger, those
four `outcome.md` files and this file.

Evidence: run 45598541744 (base e0747fc1, head cc3d4bdc; 92 test files and R CMD check): 84 files
the same, the 8 differences each as predicted, the check identical on both sides (0 errors, 3
warnings, 3 notes; 238 paths in the portable-file-names note on each side); MOFA lane run
45598538505 on the same head: `test-integration-workbook-mofa.R` 49 passed, 0 skipped;
`test-integration-record-mofa.R` 29 passed, 0 skipped. 578.1 of 2,020 job-minutes and 14.6 of 60.

Note for readers of VAL-256: on the same evening the owner changed the duplicate-resolution rule
KIT-214 implements (VAL-254: the most completely measured row first, then the mean); that change is
ticket KIT-258, with its own replay marker, and is not in this merge.

Not shown by this record: the full suite on the landed tree (next checkpoint).
