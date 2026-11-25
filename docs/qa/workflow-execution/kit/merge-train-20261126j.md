# Merge train 17 of 2026-11-26: KIT-246

**This train is STACKED on merge train 16, itself stacked on merge train 15.** It was branched from
`origin/merge-train-20261126i` at its tip `0fd29a6bd65f509483151492a916c7f6f895249a` (train 16's declaration commit,
after its merge 5ff6cfcf of KIT-261), which is stacked on train 15's tip 069c2832 (`merge-train-20261126h`). Both
runs were out or unlanded when this train was built; staging was d7cfe160 (trains 13 and 14 landed; ledger rows through
VAL-301). **Train 17 lands only if trains 15 and 16 both land with their code unchanged.** If either is rebuilt or its
code changes, this train is rebuilt and its check repeated.

## Train

Base: `0fd29a6bd65f509483151492a916c7f6f895249a` (train 16's tip, fetched before building).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-246 | `origin/kit-246b` d9b17540 | f4ba3a42 | `VAL-<train17-1>`; INT-31 left as it is | AUDIT-KIT-246-impl-20261126 (verdict ready; C1 to C3) |

Head (bound by the run): `f4ba3a42` (full hash in the request). This file, the predictions and the `<commit>` fill
below are in the declaration commit after it, under `docs/` only; the records commit after the run likewise.

### Conflicts met and their resolutions

- **`docs/methods-decision-ledger.md`**: KIT-246's `VAL-<merge>` row against the base's appended rows (VAL-295 to
  VAL-297 and the train 13 to 16 placeholders). Kept the base's rows, then KIT-246's row as `VAL-<train17-1>`,
  **with the register's I03 replacement text** (C1: the claim limited to records whose model used per-feature
  scaling; VAL-251 named as deciding only the help's last sentence). Its `<commit>` is filled with the merge commit
  f4ba3a42 in the declaration commit (a commit cannot name itself). INT-31 is untouched (it keeps its own `<commit>`).
  The placeholder was also replaced in `kit-246/outcome.md` (lines 16, 94, 164).
- **`NEWS.md`**: a conflict was met, although NEWS was expected to merge cleanly: KIT-261's entry (base) and KIT-246's
  entry insert at the same place under `### Fixes`. Both kept, the base's first; no text changed. It is not under
  `R/`, `tests/` or `inst/`.
- **C2 (records)**: I01's one-sentence replacement applied in `kit-246/wording.md`: "The owner's decision (VAL-293)
  is the help's passage in full (its last sentence decided in VAL-251; the sentence before it is KIT-239's text,
  VAL-255)".
- Not done in the train, as instructed: O1 (the four phrases; parked for the owner; the passage merges as built) and
  C5 (the comment at `tests/testthat/test-integration-results.R:599`).
- No conflict under `R/`, `tests/` or `inst/`. `DESCRIPTION` and `NAMESPACE` unchanged. Files changed against the base
  outside `docs/`: `NEWS.md`, `R/utils_integration_gui.R`, `R/utils_integration_record.R`,
  `tests/testthat/test-integration-record.R`.

### Path lengths

Checked by script before sending: 13 paths are new since the base; none exceeds 100 bytes with the `MultiScholaR/`
prefix; the longest is 72 (`docs/qa/workflow-execution/kit/kit-246/r2-failures-base.txt`).

## Declaration (written before the run is sent)

**Purpose.** KIT-246's merge check on the merged tree (C3).

**Scope.** 11 files, base and head, R CMD check on both sides: R2's list (`kit246-r2-c74f2a95`):
`test-integration-record.R`, `-enrichment.R`, `-outputs.R`, `-input-matrix.R`, `-workbook-api.R`,
`-removed-factors.R`, `-model.R`, `-enrichment-mofa.R`, `-results.R`, `-gui.R`, `test-kit240-integration-default.R`.

Recount on this base: the ticket's grep (`grep -lE "integrationMethodsText|integrationStudyParameterLines|
integrationRecordTexts|methods_scaling" tests/testthat/test-*.R`) gives the same 12 hits at the base 0fd29a6b and at
the head; less the ticket's four excluded readers of other keys (`-enrichment-tests.R`, `-enrichment-options.R`,
`-outputs-mofa.R`, `-record-mofa.R`) that is 8 files, all in the 11; with `-gui.R` and `-results.R` (item 2 (ii))
that is 10. The scope has not grown. One difference from the ticket's list: `test-integration-workbook-api.R` no
longer matches the grep; it is kept, as R2's declared file.

**Run settings.** `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` on both sides (`test-prot-norm-core.R` is not in scope); no extra
packages (as R2); reference BLAS; 1 shard per side (as R2); no scripts (no identity script, no MOFA lane).

**Ceiling.** 280 job-minutes: KIT-246's R2 row (merge check), one shard per side (2 × 80) and check on both (2 × 60),
ticket section 11. R2 used 115.4.

**Base values.** Trains 15 and 16 move no file in this scope (train 15 declares all eight of its overlapping files
`same` at train 14's heads; train 16's scope holds none of them); so the base values are the audit's section 4.8
values, which equal those heads.

### Predictions (AUDIT-KIT-246-impl section 4.8)

| File | Base p/f/e/s | Head p/f/e/s | Verdict |
|---|---|---|---|
| test-integration-record.R | 167/0/0/0 | 200/0/0/0 | count-changed (the KIT-246 block: passed NA -> 33) |
| test-integration-enrichment.R | 349/0/0/0 | 349/0/0/0 | same |
| test-integration-outputs.R | 183/0/0/0 | 183/0/0/0 | same |
| test-integration-input-matrix.R | 84/0/0/0 | 84/0/0/0 | same |
| test-integration-workbook-api.R | 144/0/0/0 | 144/0/0/0 | same |
| test-integration-removed-factors.R | 520/0/0/3 | 520/0/0/3 | same |
| test-integration-model.R | 265/0/0/0 | 265/0/0/0 | same |
| test-integration-enrichment-mofa.R | 0/0/0/3 | 0/0/0/3 | same |
| test-integration-results.R | 258/0/0/0 | 258/0/0/0 | same |
| test-integration-gui.R | 405/0/0/0 | 405/0/0/0 | same |
| test-kit240-integration-default.R | 92/0/0/1 | 92/0/0/1 | same |
| R CMD check | | | identical findings, nothing only on head |

### Decision and attempts

- Every prediction holds: TRAIN GREEN, landing only after trains 15 and 16 land with their code unchanged.
- Any other verdict is a result: reported with its files and counts, not repaired here.
- One run; one re-send only for an instrument defect with no result judged. An HTTP 500/502 at dispatch with no run
  created: started once by hand with the README's command, recorded.

## Result

(filled after the run)
