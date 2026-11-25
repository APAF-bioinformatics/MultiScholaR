# Merge train 14 of 2026-11-26: KIT-248

**This train is STACKED on merge train 13.** It was branched from `origin/merge-train-20261126f` at its tip
`f8c871f3ce42fc72315ab97550360bc9ce132f09` (train 13's declaration commit, after its merge a32edeb1 of KIT-254).
Train 13 is itself stacked on train 12's eb6c4000 (KIT-249, KIT-236), which is stacked on train 11 (landed as
0a774e0b). Train 13's run was out when this train was built; train 12 is not yet landed (its run's one deviation is
with KIT-249's auditor). Staging was 12701142. **Train 14 lands only if trains 12 and 13 land with their code
unchanged.** If either is rebuilt or its code changes, this train is rebuilt and its check repeated.

## Train

Base: `f8c871f3ce42fc72315ab97550360bc9ce132f09` (train 13's tip, fetched before building).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-248 | `origin/kit-248` 5dd2a6c2 | d1018f6b | `VAL-300` (`VAL-301` reserved for the coordinator's stored-value row; not written here) | AUDIT-KIT-248-impl-20261126, I03 |

Head (bound by the run): `d1018f6b` (full hash in the request). This file, the predictions table and the records
commit after the run are under `docs/qa/` only.

The merge commit's one-line summary is shorter than KIT-248's NEWS entry, which is the accurate statement: a skip of
the proteomics correlation filter is recorded as a skip ('pearson-skip', skipped = TRUE) in the analysis record and the
study parameters, and a replay keeps every sample.

### Conflicts met and their resolutions

- **No conflict under `R/`, `tests/` or `inst/`.** `R/utils_workflow_analysis_record.R` (KIT-254's replay statement
  in the base, KIT-248's pearson-skip replay branch) and `R/utils_workflow_operations.R` merged automatically.
  `R/utils_workflow_dia_products.R`: `kit-248` already held staging 0a774e0b (KIT-257) with I03's authorised
  resolution of `runDIAReferenceCorrelation`; the base holds KIT-257 too, so it merged automatically and the merged
  function is I03's text exactly (the file equals `kit-248`'s).
- **`docs/methods-decision-ledger.md`.** The base held train 11's row as `VAL-281`; `kit-248` held it landed
  as VAL-281 (identical but for the landing numbering: `VAL-282` -> VAL-282 and the run id), VAL-282 and its own
  `VAL-<merge>`. Resolved: VAL-281 and VAL-282 (staging's landed rows; VAL-281 replaces the placeholder copy, so the
  KIT-257 row is not listed twice), then the base's `VAL-295`, `VAL-297`, `VAL-298`,
  `VAL-299` unchanged, then KIT-248's row as `VAL-300`; also in `kit-248/outcome.md:237`. Staging's
  later rows VAL-283 to VAL-294 are on neither side of this merge; they are the coordinator's at landing.
- **`NEWS.md`.** The base's entries (KIT-236, KIT-254) then KIT-248's, both kept.
- `DESCRIPTION` and `NAMESPACE` unchanged.

### Path lengths

Checked by script before sending: 29 paths are new since the base; none exceeds 100 bytes with the `MultiScholaR/`
prefix; the longest is 83 (`docs/qa/workflow-execution/kit/merge-train-20261126d/check-summary.txt`).

## Declaration (written before the run is sent)

**Purpose.** KIT-248's merge check on the merged tree (register I03, "Closure evidence").

**Scope.** 81 files, base and head, R CMD check on both sides: the 80 files of R2 (`kit248-r2-29072595`, run
45649121588) and `test-kit257-norm-own-files.R`.

**Run settings.** `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` and `MULTISCHOLAR_KIT105_W3_DUMP=1` on both sides;
`extra_packages` limpa and mixOmics on both sides (I03: without them E6, E10, E10b, the limpa half of E11 and
`test-kit244-limpa-pass-through.R` skip, an instrument defect; R2's request listed no extra packages, its library
held limpa 1.4.2 and mixOmics 6.36.0); reference BLAS; 8 shards per side (files ordered by recorded seconds); per-file
kill limits of 120 minutes for `test-prot-norm-core.R` and 90 for `test-workflow-saved-da.R` (not thresholds).

**Ceiling.** 1,460 job-minutes = 16 test jobs × 80 + 2 check jobs × 90. Derivation: the ticket's R2 row is
6 × 70 + 2 × 90 = 600 at 3 shards; R2 itself ran at 4 shards and took 793.4 job-minutes (recorded overrun,
`kit-248/outcome.md`), about 383 minutes of test time per side by the recorded seconds of these 81 files. At 8 shards
that is about 48 test minutes per job plus about 15 of set-up and package installation, so 80 per test job leaves a
margin; the check jobs keep the ticket's 90. Expected about 2 × (383 + 8 × 15) + 30 ≈ 1,040.

**Base values.** Per file in `merge-train-20261126g-predictions.tsv`: the declared head value of the latest train
whose scope holds the file (train 13, else train 12, else train 11's table at the checkpoint 45644467847), else the run
tables of trains 10, 9, 8. Train 13 moves no file here except via its base (`test-omics-norm-state-guards.R` 651/0);
train 12 moves `test-workflow-saved-da.R` to 755/196 (observed so in run 45667562462); train 11 gives
`test-kit257-norm-own-files.R` 133/0 and `test-prot-norm-core.R`'s four KIT-257 labels. Three files are on no train
record and are predicted `same` with their counts read from the run: `test-kit172-log2-rule.R`,
`test-module-ci-trx-qc.R`, `test-workflow-omics-owner.R` (KIT-248 R2 had 212/0, 34/1/1 error, 78/1); and
`test-workflow-partitions.R`'s base is likewise read from the run. The checkpoint and train 8 installed packages
this run does not (`sva`); this run installs limpa and mixOmics, which earlier trains 10 to 13 did not, so skip counts
in limpa or mixOmics files may read differently from those records; the verdict is base against head in this run.

### Predictions (register I03)

- `test-kit248-correlation-skip-record.R`: only on head, 206 passed, 0 failed, 0 skipped. Any skip is an instrument
  defect.
- `test-prot-norm-core.R`: 1548/306 on both sides, `message-changed` at exactly 12 labels: the provenance and record
  digests of dia_correlation_skip, dia_normalization_refused and dia_normalization_retry, in W1 (a) and W2 (base
  `4aac9383…`/`a14c1ee3…`, head `dcfdf833…`/`6844e105…`); the `files` labels of refused and retry carry KIT-257's values
  on both sides. Nothing else in the file moves.
- `test-workflow-operation-registry.R`: 1095 -> 1099 passed (the KIT-039 labels loop 354 -> 356, the KIT-027 defaults
  loop 264 -> 265, the facts loop 206 -> 207); 2 failures on both sides, each only the `integration_mofa_model`
  difference (verdict `message-changed` or `count-changed`, judged on that content).
- `test-workflow-partitions.R`: 439 -> 440 passed (the compute-role loop), `count-changed`.
- Every other file `same` at the base's counts (for example `test-workflow-saved-da.R` 755/196,
  `test-kit257-norm-own-files.R` 133/0, `test-kit244-limpa-pass-through.R` 204/0 with 0 skipped,
  `test-omics-norm-state-guards.R` 651/0).
- R CMD check: findings identical on both sides, nothing only on head; the portable-file-names NOTE the same count on
  both sides.

### Decision and attempts

- Every prediction holds: TRAIN GREEN, landing only after trains 12 and 13 land with their code unchanged.
- Any other verdict is a result: reported with its files and counts, not repaired here.
- One run; one re-send only for an instrument defect with no result judged. An HTTP 500/502 at dispatch with no run
  created: started once by hand with the README's command, recorded.

## Result

Run 45676428221 (request `merge-train-20261126g-d1018f6b`, dispatched normally; `meta.txt` names the request file,
package_commit d1018f6b, base_commit f8c871f3; its request_commit reads 3b2fc169, another agent's later request commit
on `run-tests`, which contains this request file unchanged), sent once; conclusion success; 809.1 job-minutes of the
1,460 ceiling (tests 778.1, check 30.3, plan and results 0.8). Library: limpa 1.4.2, mixOmics 6.36.0 installed. Tables
in `merge-train-20261126g/` (`compare.tsv`, `tests.tsv`, `tests-base.tsv`, `check-compare.txt`, `check-summary.txt`,
`meta.txt`).

Observed, 81 files: 77 `same`, 4 not `same`:

| file | base p/f/e/s | head p/f/e/s | verdict | prediction |
|---|---|---|---|---|
| test-kit248-correlation-skip-record.R | absent | 206/0/0/0 | only-on-head | held (0 skipped) |
| test-prot-norm-core.R | 1548/306/0/0 | 1548/306/0/0 | message-changed | held: exactly the 12 labels (below) |
| test-workflow-operation-registry.R | 1095/2/0/0 | 1099/2/0/0 | count-changed | held: the 2 failures identical on both sides, each the `integration_mofa_model` difference |
| test-workflow-partitions.R | 439/0/0/2 | 440/0/0/2 | count-changed | held (base read from the run: 439, as R2) |

`test-prot-norm-core.R`, parsed from both failure files (keyed by test and leaf; 282 keys on each side, none only on
one side): 12 changed messages, exactly the provenance and record leaves of dia_correlation_skip,
dia_normalization_refused and dia_normalization_retry in W1 (a) and W2 (for example dia_correlation_skip provenance:
actual `4aac9383…` on base, `dcfdf833…` on head). Nothing else in the file moved.

Every other file `same` at the base counts on record (0 differences, script over `compare.tsv` against
`merge-train-20261126g-predictions.tsv`), including `test-workflow-saved-da.R` 755/196, `test-kit257-norm-own-files.R`
133/0, `test-kit244-limpa-pass-through.R` 204/0 with 0 skipped, `test-omics-norm-state-guards.R` 651/0. Read from the
run: `test-kit172-log2-rule.R` 212/0, `test-module-ci-trx-qc.R` 38/0/1 error (KIT-248 R2 had 34/1/1 without these
packages), `test-workflow-omics-owner.R` 78/1, all `same`.

R CMD check: 0 ERRORs, 3 WARNINGs, 3 NOTEs on both sides; all 6 findings identical, none only on head or only on base.

Verdict: every prediction held. TRAIN GREEN, landing only after trains 12 and 13 land with their code unchanged.

Not checked: trains 12 and 13 at landing; any file outside the 81; the full suite (KIT-248's A4 at the next
checkpoint).

## Merge record (coordinator)

Landed on `v0.9-staging` on 2026-11-26 by a merge of `merge-train-20261126g` (tip 0c8866ec) into
staging 988512e6. The stacking condition held: trains 12 and 13 landed with their code unchanged
(KIT-249 = VAL-295, KIT-236 = VAL-297, KIT-254 = VAL-298 and VAL-299), and staging differs from
the commit this train was built on (f8c871f3) only under `docs/`. The only conflict at landing was
the ledger: staging's rows kept through VAL-299, then this train's row. Ledger rows: KIT-248 =
VAL-300; the coordinator's stored-value row for the six `provenance` and `record` values of
`tests/testthat/test-prot-norm-core.R` (already on the re-record list since KIT-187, VAL-104; not
edited) = VAL-301. In the same commit the line numbers in VAL-282 are corrected (768 and 802; the
block moved by one line with KIT-244's re-pin).

Evidence: run 45676428221 (base f8c871f3, head d1018f6b; 81 test files and R CMD check, limpa and
mixOmics installed, both environment switches set on both sides): 77 files the same; the new test
file only on head and passing (206/0, none skipped); `test-prot-norm-core.R` 1548/306 on both
sides with exactly the twelve declared messages changed; the registry and partitions files with
their declared count changes; the check identical on both sides (0 errors, 3 warnings, 3 notes).
809.1 of 1,460 job-minutes.

With this landing KIT-244 and KIT-248 are both on staging: on DIA with limpa, "press again, then
Skip" no longer gives a live result that differs from replay.

Not shown by this record: the full suite on the landed tree (next checkpoint).
