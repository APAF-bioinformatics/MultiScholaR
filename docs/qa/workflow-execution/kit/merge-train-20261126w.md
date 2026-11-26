# Merge train 30 of 2026-11-26: KIT-258 texts (VAL-306)

**This train is NOT stacked on other open trains.** It was branched from `origin/v0.9-staging` at
`9171349894d7212b37c29281cb11447e22c8d59e` (fetched before building; staging's newest at the time; no train carrying
KIT-278 has landed on it). Base counts are therefore those of the register's table (R2 and R3 of the branch, and the later
committed tables), and the train's own base run is the reference where they differ.

## Train

Base: `9171349894d7212b37c29281cb11447e22c8d59e` (staging).

| Order | Ticket | Branch head | Merge commit | Ledger rows | Audit |
|---|---|---|---|---|---|
| 1 | KIT-258 texts (VAL-306) | `origin/kit-258b` 97e89727 (its base 75acd421) | bc62602f | `VAL-<train30-1>` (the texts) and `VAL-<train30-2>` (stored references), numbered by the coordinator at landing | AUDIT-KIT-258-impl-20261126, sections "Texts (VAL-306)" and "Closure check (fifth text)": READY FOR A MERGE TRAIN |

Head (bound by the run): the last code-bearing commit of this branch, which is the merge commit `bc62602f49a730d954d4dbd48821711f20d9d5b3`
(the record nits C1 and C2 and the two ledger rows are in it; documents only). This file is in the declaration commit after it,
under `docs/` only.

### Conflicts met and their resolutions

- **`docs/methods-decision-ledger.md`** (the only conflict): every staging row kept (up to VAL-354); directly after VAL-354
  the branch's row as `VAL-<train30-1>` (its id cell was `VAL-<merge>`), then the stored-reference row as `VAL-<train30-2>`,
  with the exact text of the register's closure check (section (d), "Final stored-reference row", its id cell `VAL-<merge-b>`
  written `VAL-<train30-2>`; the register's wrapped lines joined with single spaces).
- `R/func_metab_s4_qc_methods.R` auto-merged (staging changed its correlation histogram lines since 75acd421).
- No conflict under `R/`, `tests/`, `inst/`, `man/`, `DESCRIPTION` or `NAMESPACE`.

### Merged code against the base

`git diff --stat origin/v0.9-staging bc62602f -- R tests inst man DESCRIPTION NAMESPACE`: 12 files, 177 insertions, 20
deletions, exactly the branch's own change against its base 75acd421 (the same 12 files, 177/20; the changed lines
(`+`/`-` lines of the two diffs) are identical, compared by script):
- text strings: `R/func_lipid_qc_support_helpers.R` (6), `R/func_lipid_s4_duplicate_helpers.R` (4),
  `R/func_metab_s4_qc_methods.R` (6), `R/func_omics_qc_shared.R` (6), `R/mod_omic_qc_tabs.R` (3),
  `R/utils_omic_gui_config.R` (6): the help, the result text, the console messages, two roxygen titles, one comment;
- the step description: `R/utils_workflow_omics_declarations.R` (2, the literal at `:1853`);
- a comment: `R/utils_workflow_lipid_declarations.R` (5);
- the two `\title` lines: `man/resolveDuplicateFeaturesByIntensity.Rd`, `man/resolveLipidDuplicateFeaturesByIntensity.Rd` (2 each);
- the new test `tests/testthat/test-kit258b-duplicate-texts.R` (153) and the one pin edit
  `tests/testthat/test-lipid-01-qc-filtering-helpers.R` (2).
Nothing under `inst/`, no `DESCRIPTION`, no `NAMESPACE`; no other test file.

### Record fixes in the merge commit (documents only)

- C1: `kit-258/outcome.md`, "Not checked" paragraph: "(both done since: see the Addition below)" appended.
- C2: `kit-258/texts-step0.md` section 5: "the file reads 40 passed, 6 failed (the six of block 5)".

## Declaration (written before the run is sent)

**Purpose.** The merge check of KIT-258's in-app texts on the merged tree: base against head on the register's scope, so
that the declared message changes are confirmed on the landing base. The full suite is not repeated.

**Scope recount (on the merged tree).** The register's 35 files: R2's 34 (request `kit258b-r2-ac12c7e0`, quarantined
files out) plus `test-workflow-omics-binding-steps.R` (the step description). All 35 exist on the merged tree; the new file
is the only one absent on the base. Recount by `git grep -lE` over `tests` (the function names, `duplicate_resolution`,
`omicQcDuplicatesText`, `duplicates_help`, `qc-duplicates-resolution`, the two QC duplicates modules, "Resolved duplicate",
"highest average intensity", `qc_duplicates`; `data/benchmark` excluded by path): 40 hits, of which 26 are in the 35; the
other 14 are not test files to run (three helper files, the four gui-entry JSON data files and `tests/testdata/mci/manifest.json`,
`tests/verify_fixes.R`) or are quarantined (`skip("quarantined: Q-56/57/60/66/67`: `test-lipid-01b`, `01c`, `01f`,
`test-metab-01o`, `01p`); the 9 files of the 35 that the grep does not name are the GUI-entry lanes, the `01h`/`01t` S4
characterisations, the globals guard, the QC core, the harness, the residency test and the binding-steps test (the
register's step 0 names them by their readers).

**Settings.** As the branch's R2 request (`kit258b-r2-ac12c7e0`): `env` `{}`, `blas` reference, `check` true, no extra
packages, `file_timeout_minutes` 45, no `file_timeouts`, **2 shards per side**. Base = `9171349894d7212b37c29281cb11447e22c8d59e`,
head = `bc62602f49a730d954d4dbd48821711f20d9d5b3`; one run.

**Ceiling.** 220 job-minutes = 14.67 core-hours (core-hours = job-minutes x 4 / 60), R2's ceiling. Expected about 145
job-minutes (R2's 141.07 plus about 2 to 3 per side for the binding-steps test) = about 9.7 core-hours.

### Predictions

Base counts from the register (R2, R3 and later committed tables: `merge-train-20261126p`, `m`, `h`, `kit-264/r2`);
INFERRED for files no later run covered; the train's own base run is the reference. Counts are passed/failed
(skipped in brackets where the table has them).

| File | Base | Head | Verdict |
|---|---|---|---|
| test-kit258b-duplicate-texts.R | absent | 46/0 | only-on-head |
| test-workflow-omics-binding-steps.R | 57/2 | 57/2 | message-changed (exactly the two `steps[[4]]$description` lines of the KIT-094 snapshot, lipidomics and metabolomics) |
| test-omic-qc-module.R | 156/4 | 156/4 | message-changed (`tabs` and `subtabs` screen digests; `ids` and `id_counts` identical) |
| test-gui-entry-lipidomics.R | 19/16 | 19/16 | message-changed (the 15 steps from qc_duplicates on) |
| test-gui-entry-metabolomics.R | 57/48 | 57/48 | message-changed (45 steps: 15 in each of three lanes) |
| test-lipid-01-qc-filtering-helpers.R | 58/0 (2) | 58/0 (2) | same |
| test-omics-module-globals.R | 13/1 | 13/1 | same (14/0 on both sides if the base run shows KIT-278's change; none on this base) |
| test-omics-norm-state-guards.R | 651/0 | as base | same |
| test-workflow-omics-worker.R | 289/0 | as base | same |
| test-workflow-partitions.R | 440/0 (2) | as base | same |
| test-workflow-omics-residency.R | 97/0 | as base | same |
| test-kit258-duplicate-completeness.R | 162/0 | as base | same |
| test-kit214-sample-columns.R | 176/0 | as base | same |
| test-lipid-01g-qc-duplicates-module-characterization.R | 9/0 (5) | as base | same |
| test-lipid-01h-qc-s4-module-characterization.R | 3/0 | as base | same |
| test-lipid-01j-qc-orchestrator-shared.R | 3/0 (3) | as base | same |
| test-lipid-03-s4-duplicate-helpers.R | 26/0 | as base | same |
| test-metab-01b-qc-filtering-helpers-characterization.R | 9/0 | as base | same |
| test-metab-01m-s4-resolve-duplicate-intensity-characterization.R | 19/0 | as base | same |
| test-metab-01t-qc-s4-module-characterization.R | 2/0 (7) | as base | same |
| test-metab-01u-qc-duplicates-module-characterization.R | 9/0 (3) | as base | same |
| test-metab-01w-qc-orchestrator-shared.R | 3/0 (3) | as base | same |
| test-metab-02u-resolve-duplicate-features-characterization.R | 17/0 | as base | same |
| test-metab-s4-qc-shared.R | 32/0 | as base | same |
| test-module-ci-lipid-qc.R | 68/0 | as base | same |
| test-module-ci-metab-qc.R | 63/0 (1) | as base | same |
| test-module-ui-wrapper-direct-shared.R | 45/0 | as base | same |
| test-omic-gui-config.R | 1911/0 | as base | same |
| test-omics-audit-carry-overs.R | 20/0 | as base | same |
| test-omics-qc-core.R | 115/0 | as base | same |
| test-workflow-kit-harness.R | 734/1 | as base | same |
| test-workflow-lipid-declarations.R | 259/0 | as base | same |
| test-workflow-metab-declarations.R | 259/0 | as base | same |
| test-workflow-omics-declarations.R | 43/25 | as base | same |
| test-workflow-operation-registry.R | 1107/2 | as base | same |
| R CMD check | 0 E / 3 W / 3 N | identical | identical findings, nothing only on head |

The 35 files: 1 `only-on-head`, 4 `message-changed` (binding-steps, QC module, the two GUI-entry lanes), 30 `same` (the register's
"other 28" plus `test-lipid-01-qc-filtering-helpers.R` and `test-omics-module-globals.R`).

### Verification of the masked comparison (planned)

Script: Python `-I`, kept outside the downloaded folder; input the run's `failures.txt` (head) and `failures-base.txt`
(base), records split at `file:`/`test:`/`kind:`/`message:`; the output is committed. Step 1: for every file, the records
that differ between base and head after the run's own path masking are listed by file and test. Step 2, the three
GUI-entry files (lipidomics, metabolomics with its three lanes, and the totals): every record outside the 15 steps from
`qc_duplicates` on must be identical on both sides; for each of those 15 steps per lane, the base record and the head record
are compared line by line (`difflib`): the only lines that differ must contain `Strategy: Keep feature with highest
average intensity` (base) or `Rule:` followed by the new sentence (head), and the record is then equal to its base record
after those lines are deleted; any other residue is a missed prediction. The 4 other failing records per file stay identical
(60 differing records in all, 15 x 4 lanes). Step 3, `test-omic-qc-module.R`: the records of `tabs` and `subtabs` differ
by the digest value only (`actual`); the `ids` and `id_counts` records are identical. Step 4, `test-workflow-omics-binding-steps.R`:
the two records differ only in the lines naming `steps[[4]]$description` (the old words against the new sentence).

### Decision and attempts

- Every prediction holds: TRAIN GREEN for KIT-258 texts.
- A missed prediction is a result: recorded, handed back, not repaired and not re-sent. A mechanics failure with no
  result: handed back with the evidence; no re-send without the coordinator's word. An HTTP 500/502 at dispatch with no
  run created: started once by hand with `gh workflow run run-tests.yml --ref run-tests -f request=<name>`.

## What the coordinator resolves at landing

- The two ledger numbers: `VAL-<train30-1>` (the texts) and `VAL-<train30-2>` (stored references): the id cells in
  `docs/methods-decision-ledger.md`, and the references to them in the two rows; `<date>` and the run id where the rows have them.
- The status cell of VAL-337 and of VAL-306 (the texts, AUDIT-KIT-258-impl-20261126 section 6 (T08)); the sentence of
  the branch's row about the step description (the coordinator's routing of T02: ruled inside VAL-306's intent, applied).
- Placeholders in the branch: see the list at the end of the run's results section.

## Result

Run 45737029722 (request `train30-kit258-texts-bc62602f`, request commit 3d245a5b), dispatched by the dispatcher without error
(no hand dispatch); green (infrastructure). Head bc62602f vs base 91713498. Tables: `merge-train-20261126w-compare.tsv`,
`-tests.tsv`, `-tests-base.tsv`, `-check-compare.txt`, `-meta.txt`, `-masked-comparison.txt` (same directory).

**Job-minutes 164.63 = 10.98 core-hours** (job-minutes x 4 / 60) against the ceiling of 220 (14.67 core-hours); expected about
145 (so over the expectation, under the ceiling): plan 0.40, tests head 42.15 and 30.10, tests base 43.37 and 26.18, check head
9.38, check base 12.48, results 0.57 (job start to end from the jobs API; skipped scripts job 0).

| File | Predicted | Observed base | Observed head | Verdict |
|---|---|---|---|---|
| test-kit258b-duplicate-texts.R | only-on-head, 46/0 | absent | 46/0 | only-on-head: HELD |
| test-workflow-omics-binding-steps.R | 57/2 both, message-changed (two description lines) | 57/2 | 57/2 | message-changed: HELD |
| test-omic-qc-module.R | 156/4 both, message-changed | 156/4 | 156/4 | message-changed: HELD |
| test-gui-entry-lipidomics.R | 19/16 both, message-changed | 19/16 | 19/16 | message-changed: HELD |
| test-gui-entry-metabolomics.R | 57/48 both, message-changed | 57/48 | 57/48 | message-changed: HELD |
| test-lipid-01-qc-filtering-helpers.R | 58/0 (2) same | 58/0 (2) | 58/0 (2) | same: HELD |
| test-omics-module-globals.R | 13/1 same | 13/1 | 13/1 | same: HELD (no KIT-278 on this base) |
| the other 28 files | same, counts as the table above | as predicted | as predicted | same: HELD (every count equals its prediction) |
| R CMD check | identical (0 E / 3 W / 3 N) | 0 E / 3 W / 3 N | 0 E / 3 W / 3 N | identical: HELD (6 findings on both; only on head: none; only on base: none) |

Totals: 35 files; 30 `same`, 4 `message-changed`, 1 `only-on-head`; no `new-failure`, no `count-changed`, no `fixed`. No
absolute-count miss.

**Masked comparison** (script: Python `-I`, kept outside the downloaded folder; output committed in
`merge-train-20261126w-masked-comparison.txt`). 99 failure records on each side; 36 identical after temporary-path masking;
63 differ, all in the four message-changed files: 45 (metabolomics GUI-entry), 15 (lipidomics GUI-entry), 2 (QC module), 1
(binding-steps). GUI-entry lanes: in each of the 60 differing records the head holds one additional chunk, the `lines(actual$outputs$<omic>_workflow-qc-duplicates-resolution_results)[1:6]` block, whose only
changed lines are `- "Rule: the row measured in the most design samples (a missing value or 0 is not measured), then the highest
average of its measured values"` against `+ "Strategy: Keep feature with highest average intensity"`; with that chunk (header, context
lines, the two changed lines and its blank line) removed, the head record **equals the base record exactly: 60 of 60, no residue**
(the 4 other failing records, one in lipidomics lane's and three in metabolomics, are identical on both sides; the 45 metabolomics
records are the 15 steps in each of three lanes). Binding-steps: the only differing lines are the two `steps[[4]]$description`
pairs (new sentence against the snapshot's "Resolved duplicate lipid/metabolite features by keeping highest intensity"). QC
module (checked by hand and by digest masking): `ids` and `id_counts` identical on both sides; `tabs` differs at `lipid_qc` and
`metab_qc` digests only, the expected values unchanged; `subtabs` differs at base in `lipid_intensity` and `metab_intensity`
(KIT-217) and on head additionally at `lipid_duplicates` and `metab_duplicates`, the expected values unchanged. Everything
differs exactly as declared.

**Decision.** Every prediction held. TRAIN GREEN for KIT-258 texts.

## Rule slips

None found. The run-tests worktree held other agents' untracked folders, untouched. One `fetch.sh`, one `gh run watch`
(backgrounded by the shell after 180 s; waited for on a local process monitor, no API polling). A few `gh run list` /
`gh run view` calls (4) were made after the push before the run appeared (the dispatcher needed about a minute). In the
declaration, the count sentence of the predictions was first garbled and is corrected in this commit.

## Placeholders for the coordinator, by file and line

- `docs/methods-decision-ledger.md` line 2527: id cell `VAL-<train30-1>`; the row's last-but-one sentence ("Not changed, found: ...
  needs a scope decision" form) is to be replaced with the coordinator's routing of T02 (applied by commit ef570448); status
  cell "implemented on `kit-258b`; not merged ..." to be set at landing.
- `docs/methods-decision-ledger.md` line 2528: id cell `VAL-<train30-2>`; status cell ("recorded; deferred to the end-of-project
  re-record") stays.
- The status cells of VAL-337 and VAL-306 (texts of the register, section (f) and T08), not edited here.
- This record: the table in "Train" (ledger numbers) and "What the coordinator resolves at landing" (the two ids).
- `docs/qa/workflow-execution/kit/kit-258/outcome.md:180`: "VAL-<merge>" in the A5 row (the branch's own text, historical).
