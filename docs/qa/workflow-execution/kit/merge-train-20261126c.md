# Merge train 10 of 2026-11-26: KIT-259

Built on the coordinator's instruction after trains 8 and 9 landed (staging 828dd105, ledger rows through VAL-272).
Not stacked: the base is staging. Same procedure as trains 8 and 9 (`merge-train-20261126a.md`,
`merge-train-20261126b.md`).

## Train

Base: `828dd1051d2acab075d6d04f887b7c9e9ccaa48e` (`origin/v0.9-staging`, fetched before building; unchanged).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-259 | `origin/kit-259` 485f722d | 4baea49f | `VAL-273` | AUDIT-KIT-259-impl-20261126, section 5 (I07), section 6 (I13) |

Head (bound by the run): `4baea49f` (full hash in the request). This file and the records commit after the run are
under `docs/qa/` only.

The merge commit's one-line summary is imprecise (it was written before reading the NEWS entries and left as written,
with no amend). KIT-259's NEWS entries are the accurate statement: the proteomics Pearson histogram shows every
within-group pair (bars over 4 pairs were dropped); proteomics with no two samples in one group says so; lipidomics
and metabolomics show correlations below 0 and count the pairs with no correlation; a proteomics QC figure that cannot
be drawn no longer leaves the normalisation unfinished, and is replaced by a note giving the reason.

### Conflicts

None: `kit-259` already contains 828dd105. Its ledger row sits after VAL-272, written `VAL-273`; the same
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

Dispatch: the request commit c277d799 was pushed once; the push answered HTTP 500 and the retry found it already on
origin, but no dispatcher run was created for it (none listed at c277d799; the next dispatcher run, another agent's at
9dfa3f8f, compares only its own push). With no run of this request in existence, it was started once by hand with the
README's command (`gh workflow run run-tests.yml --ref run-tests -f request=merge-train-20261126c-4baea49f`): run
45650328815, whose `meta.txt` names this request file, package_commit 4baea49f and base_commit 828dd105. This is the
first and only run of the request, not a re-send. The one `gh run watch` ended early (exit 1, a rate-limited read);
one status read afterwards: completed, success.

Run 45650328815: conclusion success; 487.4 job-minutes of the 820 ceiling (tests 458.2, check 28.3, plan and results
0.9). Tables in `merge-train-20261126c/` (`compare.tsv`, `tests.tsv`, `tests-base.tsv`, `check-compare.txt`,
`check-summary.txt`, `meta.txt`).

Observed, 53 files: 52 `same`, 1 `only-on-head`:

| file | base p/f/e/s | head p/f/e/s | verdict | prediction |
|---|---|---|---|---|
| test-kit259-pearson-histogram-figure-failure.R | absent | 85/0/0/0 | only-on-head | held |
| test-metab-01k-s4-plot-pearson-characterization.R | 17/0/0/0 | 17/0/0/0 | same | held |
| test-workflow-saved-exports.R | 271/0/1/4 | 271/0/1/4 | same | held |
| test-prot-norm-core.R | 1548/306/0/0 | 1548/306/0/0 | same | held (failure list identical; no failure names the pass-through body digest) |
| test-workflow-saved-da.R | 776/175/0/0 | 776/175/0/0 | same | held |
| test-gui-entry-proteomics.R | 168/25/0/0 | 168/25/0/0 | same | held |
| test-kit244-limpa-pass-through.R | 204/0/0/0 | 204/0/0/0 | same | held |

Every other file `same` at its staging count on record. The eleven not on record, read from the run (both sides):
`test-lipid-05-s4-correlation-helpers.R` 21/0, `test-metab-norm-standalone-shared.R` 21/0,
`test-omics-cyclic-loess-minimum.R` 72/0, `test-omics-norm-core.R` 62/0 (14 skipped),
`test-omics-norm-itsd-selection.R` 22/0, `test-prot-s4-qc-direct-shared.R` 72/0,
`test-prot-04d-peptide-s4-normalization-characterization.R` 112/0, `test-metab-s4-qc-shared.R` 32/0,
`test-prot-05c-norm-module-seams.R` 0/0 (4 skipped), `test-metab-01k-...` 17/0 (above), and the new file (above).

R CMD check: 0 ERRORs, 3 WARNINGs, 3 NOTEs on both sides; all 6 findings identical, none only on head or only on base.

Verdict: every prediction held. TRAIN GREEN.

Not checked: any file outside the 53; the full suite (KIT-259's A-criterion at the next checkpoint; a staging
checkpoint run by another agent was out at the same time and is not read here); KIT-259's partners (KIT-257, KIT-242,
KIT-247) merging later.

## Merge record (coordinator)

Landed on `v0.9-staging` on 2026-11-26 by a merge of `merge-train-20261126c` (tip fd6ff63c) into
staging abe2199c, which differs from the train's base 828dd105 by docs only (the full-suite
checkpoint's records and three docs-only corrections). Ledger row: KIT-259 = VAL-273 (in the
ledger, `kit-259/outcome.md` and this record), its status cell filled. No conflict at landing.

Evidence: run 45650328815 (base 828dd105, head 4baea49f; 53 test files and R CMD check, both
environment switches set on both sides): 52 files the same, the new test file only on head and
passing (85/0); the check identical on both sides (0 errors, 3 warnings, 3 notes). 487.4 of 820
job-minutes. The request's first push returned an HTTP 500 and created no run; the run was started
once by hand with the documented command (a first send, not a re-send).

Open before the release commit (the owner's, parked row 45): an inferred axis-label defect whose
repair needs a second correction, and the report text that describes a placeholder panel as
drawn. The one-line summary of merge commit 4baea49f is inaccurate; the NEWS entries describe the
change correctly.

Not shown by this record: the full suite on the landed tree (next checkpoint).
