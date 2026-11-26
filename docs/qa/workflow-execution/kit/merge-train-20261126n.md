# Merge train 21 of 2026-11-26: KIT-267

**This train is STACKED on merge train 20 (KIT-258), itself stacked on train 19 (KIT-264), which is stacked on train
18 (landed on staging as 15d86bd4 with trains 15 to 17).** It was branched from `origin/merge-train-20261126m` at its
tip `22c50184318c406f3213df929d06cc21a19c585e` (fetched before building; it carries the KIT-258 merge 77b67929 and
train 20's declaration). Staging is 176f9d3f. **Train 21 lands only if trains 19 and 20 both land with their code
unchanged.** If either is rebuilt or its code changes, this train is rebuilt and its check repeated. The runs of
trains 19 (45707095842, head 9f65a58c) and 20 (45712155865, head 77b67929) were still in progress when this train was
built; their results are not in this record and this train does not touch them.

## Train

Base: `22c50184318c406f3213df929d06cc21a19c585e` (train 20's tip; the code is that of 77b67929).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-267 | `origin/kit-267` 05c71da7 | 67e151c9 | `VAL-340` (replaces `VAL-<merge>`) | AUDIT-KIT-267-impl-20261126 (verdict ready; I09 = this plan) |

Head (bound by the run): `67e151c97caf97cbc85dc65643b32b6c6b2708fa`. This file and the predictions are in the
declaration commit after it, under `docs/` only; the records commit after the run likewise.

### Conflicts met and their resolutions

- **`docs/methods-decision-ledger.md`** (the only conflict): the base's appended rows (`VAL-336`,
  `VAL-337`) kept in order, then staging's ledger-only rows VAL-300 to VAL-313 (the three commits the trains
  lack: 82a5df02, f49b83a8, 4d588ba3) exactly as staging has them, then KIT-267's row as `VAL-340` with the
  audit's I01 text (cell 4: "implementation (KIT-267) under VAL-310 ...; runs 45704530515 and 45707020223;
  AUDIT-KIT-267-impl-20261126"; cell 5: "merged into `v0.9-staging` 2026-11-26 by merge train 21; the organism mark's
  words are a draft for the wording pass (VAL-41)"). VAL-310's status cell: "decided; implemented by KIT-267
  (VAL-340), merged by merge train 21". Nothing inherited was renumbered; the base's placeholders are NOT
  repaired here.
- **`NEWS.md`** and **`R/utils_integration_gui.R`**: merged automatically. KIT-267's NEWS line (the Enrichment tab's
  defaults) is rewrapped at about 120 characters in the merge commit, no word changed (audit I07).
- **`docs/qa/workflow-execution/kit/kit-267/outcome.md`** lines 43 to 46: the cause of R1's overrun replaced with the
  audit's I02 sentence (a records edit in the merge commit).
- No conflict under `R/`, `tests/`, `inst/`, `man/`, `DESCRIPTION` or `NAMESPACE`.

### Path lengths

Checked before sending: every path new or changed since the base is at most 72 bytes with the `MultiScholaR/` prefix
(`docs/qa/workflow-execution/kit/kit-267/r2-check-compare.txt`). Longer paths in the tree (88 to 97 bytes, under
`kit-014` and `analysis-record-series`) are older and the same on both sides.

## Declaration (written before the run is sent)

**Purpose.** KIT-267's merge check on the merged tree.

**Scope recount.** The ticket's grep (`git grep -lE "integrationEnrichmentParameterSpec|integrationSettings|
mod_integration_enrichment_ui|integrationParameterProvisional|mark_provisional|mark_decided|integrationEnrichmentLabel|
integrationGuiTexts|integrationRecordTexts|settled_by" -- 'tests/testthat/test-*.R'`) gives 21 files on the base
22c50184 and 22 on the head (the new `test-kit267-enrichment-defaults.R`). Scope = those **22 test files**, the same
list as R2's request and the audit's I09 expected (22): the 20 `test-integration-*` files, `test-kit240-integration-
default.R`, `test-kit267-enrichment-defaults.R`. No difference from 22.

**Run settings.** No env switches; no extra packages (the `*-mofa.R` files skip as in R2); reference BLAS; ubuntu-22.04;
4 shards per side; file timeout 45 minutes; check on both sides; no scripts (neither
`kit-240/default-identity.R`, retired by KIT-267, nor kit-138's).

**Ceiling.** 210 job-minutes = 14.0 core-hours: KIT-267's R2 row. Expected about 190 job-minutes (about 12.7
core-hours).

### Predictions

Counts are passed/failed/errors/skipped.

| File | Base | Head | Verdict |
|---|---|---|---|
| test-kit267-enrichment-defaults.R | absent | 3/0/0/0 | only-on-head |
| test-kit240-integration-default.R | 92/0/0/1 | 96/0/0/1 | count-changed, no failure on either side |
| test-integration-record.R | 200/0/0/0 | 200/0/0/0 | same |
| test-integration-enrichment.R | 349 | 349 | same |
| test-integration-enrichment-options.R | 271 | 271 | same |
| test-integration-enrichment-tests.R | 181/0/0/2 | 181/0/0/2 | same |
| test-integration-removed-factors.R | 520/0/0/3 | 520/0/0/3 | same |
| test-integration-gui.R | 405 | 405 | same |
| test-integration-model.R | 265 | 265 | same |
| test-integration-results.R | 258 | 258 | same |
| test-integration-outputs.R | 183 | 183 | same |
| test-integration-partitioned-alignment.R | 253 | 253 | same |
| test-integration-workbook-api.R | 144 | 144 | same |
| test-integration-input-matrix.R | 84 | 84 | same |
| test-integration-sample-match.R | 76 | 76 | same |
| test-integration-workbook.R | 43 | 43 | same |
| test-integration-enrichment-mofa.R | 0/0/0/3 | 0/0/0/3 | same |
| test-integration-model-tab-mofa.R | 0/0/0/2 | 0/0/0/2 | same |
| test-integration-outputs-mofa.R, -record-mofa.R, -results-mofa.R | 0/0/0/1 each | same | same |
| test-integration-workbook-mofa.R | 0/0/0/3 | 0/0/0/3 | same |
| R CMD check | | | identical findings, none only on head |

(The `kit240` change: T6 13 to 21 passes; T10's old title, 12 passed, absent on head and its new title 8 passed. If the
compare tool names the renamed T10 differently, that is the tool's naming: judge by the counts and the empty failure
list.)

### Decision and attempts

- Every prediction holds: TRAIN GREEN, landing only after trains 19 and 20 land with their code unchanged.
- Any other verdict is a result: recorded with its files and counts, not repaired or re-sent.
- One run. One re-send only for an instrument defect with no result judged. An HTTP 500/502 at dispatch with no run
  created: started once by hand with the README's command, recorded.

## What the coordinator resolves at landing

- The ledger: `VAL-340` takes its number after the rows that landed before it (staging ends at VAL-329; trains
  19 and 20 carry `VAL-336` and `VAL-337`); `<t>` and the VAL-310 status cell's number fill. The
  inherited placeholders are reconciled as trains 19 and 20 leave them.
- `NEWS.md`: KIT-267's line sits next to KIT-246's block; order at landing.
- Audit I05 (partner edits for the retired `kit-240/default-identity.R`), I08 (optional status notes on the KIT-141
  rows) and the ticket closure (the R1 breach of its 25-minute row, calibration for later R1 rows): the coordinator's.

## Result

Run 45714416692 (request `merge-train-20261126n-67e151c9`, request commit 8e56684f), dispatched without error; green
(infrastructure). Head 67e151c9 vs base 22c50184. Tables: `merge-train-20261126n-compare.tsv`,
`merge-train-20261126n-check-compare.txt` (same directory). The first `gh run watch` ended on an API 403 (rate limit),
not on a run result; the run was read after it completed.

**Job-minutes 199.4 = 13.30 core-hours** against the ceiling of 210 job-minutes (14.0 core-hours): tests head
30.1/27.0/16.9/14.3, tests base 30.2/22.5/16.5/18.1, check head 10.4, check base 12.8, plan 0.3, results 0.5.

| File | Predicted | Observed base | Observed head | Verdict |
|---|---|---|---|---|
| test-kit267-enrichment-defaults.R | only-on-head, 3/0/0/0 | absent | 3/0/0/0 | only-on-head: HELD |
| test-kit240-integration-default.R | count-changed, 92/0/0/1 to 96/0/0/1 | 92/0/0/1 | 96/0/0/1 | count-changed (T10's old title gone, 12 passed): HELD |
| test-integration-record.R | same, 200/0/0/0 | 200/0/0/0 | 200/0/0/0 | same: HELD |
| the other 19 files (counts as predicted) | same | as predicted | as predicted | same: HELD |
| R CMD check | identical findings, none only on head | 0 E / 3 W / 3 N | 0 E / 3 W / 3 N | identical (6 findings), none only on head or base: HELD |

Totals: 22 files, 20 `same`, 1 `count-changed`, 1 `only-on-head`; no failure, error or `new-failure` on either side.

**Decision.** Every prediction held. KIT-267 introduces no regression: its new file passes 3/0 on head, KIT-240's moved
tests pass (96/0), every other file and R CMD check are identical between the train-20 tip and the merge head. TRAIN
GREEN, landing only after trains 19 and 20 land with their code unchanged (their own runs 45707095842 and 45712155865
were not finished when this train was built and are not judged here).

## Merge record (coordinator)

Landed on v0.9-staging on 2026-11-26 by a no-fast-forward merge of this branch at 7a47375e onto 0ebce762 (merge trains 19, 20 and 22 landed).

- Verdict: green. Run 45714416692 (22 files and R CMD check; 199.4 of 210 job-minutes): the new file 3/0 only on head; the kit240 file 92 to 96 count-changed with no failure; 20 files same; check identical.
- Ledger at landing: staging's rows kept through VAL-339; KIT-267 = VAL-340; VAL-310's status cell updated.
- docs/qa/workflow-execution/kit/kit-240/default-identity.R is retired from run requests from this landing on (coordinator's ruling on AUDIT-KIT-267-plan P01); the tickets and plan texts that still name it are corrected by the coordinator (AUDIT-KIT-267-impl I05).
- Recorded at closure: the ticket's first run cost 45.18 job-minutes against its step row of 25 (inside the ticket's total).
