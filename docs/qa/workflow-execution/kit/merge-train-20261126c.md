# Merge train 10 of 2026-11-26: KIT-259

Built on the coordinator's instruction after trains 8 and 9 landed (staging 828dd105, ledger rows through VAL-272).
Not stacked: the base is staging. Same procedure as trains 8 and 9 (`merge-train-20261126a.md`,
`merge-train-20261126b.md`).

## Train

Base: `828dd1051d2acab075d6d04f887b7c9e9ccaa48e` (`origin/v0.9-staging`, fetched before building; unchanged).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-259 | `origin/kit-259` 485f722d | 4baea49f | `VAL-<train10-1>` | AUDIT-KIT-259-impl-20261126, section 5 (I07), section 6 (I13) |

Head (bound by the run): `4baea49f` (full hash in the request). This file and the records commit after the run are
under `docs/qa/` only.

The merge commit's one-line summary is imprecise (it was written before reading the NEWS entries and left as written,
with no amend). KIT-259's NEWS entries are the accurate statement: the proteomics Pearson histogram shows every
within-group pair (bars over 4 pairs were dropped); proteomics with no two samples in one group says so; lipidomics
and metabolomics show correlations below 0 and count the pairs with no correlation; a proteomics QC figure that cannot
be drawn no longer leaves the normalisation unfinished, and is replaced by a note giving the reason.

### Conflicts

None: `kit-259` already contains 828dd105. Its ledger row sits after VAL-272, written `VAL-<train10-1>`; the same
placeholder replaced in `kit-259/outcome.md` (acceptance table). NEWS: KIT-259's entries under the one unreleased
heading. `DESCRIPTION` and `NAMESPACE` unchanged.

### Path lengths

Checked by script before sending: 19 paths are new since the base; none exceeds 100 bytes with the `MultiScholaR/`
prefix; the longest is 88 (`docs/qa/workflow-execution/kit/kit-259/r2c-45641244855/compare-failures.txt`).

## Declaration (written before the run is sent)

**Purpose.** KIT-259's merge check on the merged tree (register section 5; condition 2).

**Scope.** 53 files, base and head, R CMD check on both sides:
- R2's 43 files (`kit259-r2-1d6737c1`) plus `test-kit243-summary-current-outputs.R`: the ticket's two greps and five
  named files recounted at 828dd105, 44 files;
- `test-kit244-limpa-pass-through.R`;
- the eight dependents R2 never ran: `test-prot-design-log2.R`, `test-workflow-gui-memory-lifetimes.R`,
  `test-module-ci-prot-dia-limpa.R`, `test-workflow-analysis-record.R`, `test-workflow-integration-input.R`,
  `test-workflow-saved-da.R`, `test-prot-05c-norm-module-seams.R`, `test-workflow-server-shared.R`.

**Run settings.** `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` and `MULTISCHOLAR_KIT105_W3_DUMP=1` on both sides (as train 8,
for `test-prot-norm-core.R` and the KIT-244 file); no extra packages (KIT-259's own runs used none; train 8 installed
`sva` for KIT-201, train 9 did not, and train 9 read train 8's `test-gui-entry-proteomics.R` count unchanged); reference
BLAS; 6 shards per side; per-file kill limits of 120 minutes for `test-prot-norm-core.R` and the three GUI-entry files
(not thresholds).

**Ceiling.** 820 job-minutes: KIT-259's merge-check (R2) ceiling, ticket section 11 (8 × 80 + 2 × 90).

**Base values.** Staging's own counts from the latest train run on the file (`merge-train-20261126b/tests.tsv` where
the file is there, otherwise `merge-train-20261126a/tests.tsv`; trains 8 and 9 landed with only ledger edits after
their runs). Eleven files are on neither record: they are predicted `same` and their counts are read from the run
(`test-lipid-05-s4-correlation-helpers.R`, `test-metab-01k-s4-plot-pearson-characterization.R`,
`test-metab-norm-standalone-shared.R`, `test-omics-cyclic-loess-minimum.R`, `test-omics-norm-core.R`,
`test-omics-norm-itsd-selection.R`, `test-prot-s4-qc-direct-shared.R`,
`test-prot-04d-peptide-s4-normalization-characterization.R`, `test-metab-s4-qc-shared.R`,
`test-prot-05c-norm-module-seams.R`, and the new file's base side, absent).

### Predictions

- `test-kit259-pearson-histogram-figure-failure.R`: only on head, 85 passed, 0 failed, 0 skipped.
- `test-metab-01k-s4-plot-pearson-characterization.R`: `same`, 17/0 on both sides. At 828dd105 the base runs its own
  (old) test file on the old code, which passes; the head runs the edited test (`:215`) on the new code. The 16/1 in
  the correction run 45641244855 was the edited test on the old code (its base was the tests-only commit 545c4557),
  which this train's base does not contain.
- `test-workflow-saved-exports.R`: `same`, 271 passed, 0 failed, 1 error (lobstr), 4 skipped on both sides (train 8's
  count). For the same reason: the 261/10/1 base of run 45641244855 was the edited expectation (`:730-731`) on the old
  code; here the base runs staging's own expectation, which holds on the old code.
- `test-prot-norm-core.R`: `same`, 1548/306, the failure list identical, the body-digest test passing.
- `test-workflow-saved-da.R`: `same`, 776/175, messages identical.
- Every other file `same`, at the counts above where on record (for example `test-gui-entry-proteomics.R` 168/25,
  `test-kit244-limpa-pass-through.R` 204/0, `test-workflow-operation-registry.R` 1095/2,
  `test-workflow-omics-binding-steps.R` 58/1, `test-workflow-gui-diagnostic.R` 65/0/3 errors,
  `test-workflow-memory-mechanism-member.R` 13 skipped, the two quarantined files 0/0).
- R CMD check: findings identical on both sides, nothing only on head; the portable-file-names NOTE the same count on
  both sides.

### Decision and attempts

- Every prediction holds: TRAIN GREEN.
- Any other verdict is a result (KIT-259 stop rule 3): reported with its files and counts, not repaired here.
- One run; one re-send only for an instrument defect with no result judged.

## Result

(filled after the run)
