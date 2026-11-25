# Merge train 12 of 2026-11-26: KIT-249, KIT-236

**This train is STACKED on merge train 11.** It was branched from `origin/merge-train-20261126d` at its tip
`a50f93025c3aba59a2edbb41ee88431b7053cd10` (train 11's declaration commit; train 11's code head 1292afc4: KIT-257,
itself stacked on train 10's tip 53e99554). Train 11's run was out when this train was built. Staging was 3a6794b0
(train 10 landed: KIT-259 = VAL-273; ledger rows VAL-274 to VAL-280, docs only). **Train 12 lands only if train 11
lands with its code unchanged.** If train 11's check fails, or its code is rebuilt or changed, this train is rebuilt
and its check repeated.

## Train

Base: `a50f93025c3aba59a2edbb41ee88431b7053cd10` (train 11's tip, fetched before building).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-249 | `origin/kit-249` 2288d944 | 758d23ee | `VAL-295` (`VAL-296` reserved for the coordinator's correction row to VAL-237; not written here) | AUDIT-KIT-249-impl-20261125, closure check (2026-11-26) |
| 2 | KIT-236 | `origin/kit-236` 8cd81769 | 55413b59 | `VAL-297` | AUDIT-KIT-236-impl-20261125, closure check (2026-11-26) |

Head (bound by the run): `55413b59` (full hash in the request). This file, the predictions table and the records
commit after the run are under `docs/qa/` only.

### Conflicts met and their resolutions

- **KIT-249.** `kit-249` carries staging 3a6794b0, which the base (built before train 10 landed) does not. So the
  merge also brings staging's docs since train 10's tip: the ledger rows VAL-273 to VAL-280 and train 10's numbered
  records (`kit-259/outcome.md`, `merge-train-20261126c.md`), merged automatically except the ledger.
  - `docs/methods-decision-ledger.md`: the base held train 10's row as `VAL-273` and train 11's
    `VAL-281`; the branch held the same KIT-259 row landed as VAL-273 (status cell completed at landing),
    VAL-274 to VAL-280 and KIT-249's `VAL-<merge>`. Resolved: staging's landed rows VAL-273 to VAL-280 (VAL-273
    replaces the identical-but-for-status `VAL-273` row, so the row is not duplicated), then train 11's
    `VAL-281` unchanged, then KIT-249's row as `VAL-295`. This differs from "keep both" only in not
    duplicating the KIT-259 row; it is the ledger as it will read once train 11 lands.
  - `NEWS.md`: KIT-257's entry (base) then KIT-249's, both kept. KIT-249's entry already carries the closure's W2a,
    W2b and W3 corrections (I05) on the branch.
  - `R/func_prot_da_model_stats.R` and the new test merged automatically.
- **KIT-236.** `docs/methods-decision-ledger.md` (KIT-236 carried staging c4505c7f only: its row `VAL-<merge>`
  against the base's rows through `VAL-295`; kept the base's, then KIT-236's as `VAL-297`; also in
  `kit-236/outcome.md`, commit table) and `NEWS.md` (both kept, the base's first). The six report templates under
  `inst/reports/`, the KIT-227 and KIT-235 test files (KIT-236's declared edits) and the new test merged
  automatically. No conflict under `R/`, `tests/` or `inst/`.
- `DESCRIPTION` and `NAMESPACE` unchanged by the train.

### Path lengths

Checked by script before sending: 234 paths are new since the base (staging's train-10 records included); none exceeds
100 bytes with the `MultiScholaR/` prefix; the longest is 99
(`docs/qa/workflow-execution/kit/kit-236/g0-45577606908/g0a/met-src-study_parameters.txt`).

## Declaration (written before the run is sent)

**Purpose.** Both members' merge checks on the merged tree, once, on the train.

**Scope.** 103 files, deduplicated, base and head, R CMD check on both sides, and KIT-236's render job:
- KIT-249: 79 files: the 77 of its correction run (`kit249-r3-144d2c5d`: R2's 70 with the new file and the 7
  named) plus `test-kit201-treat-moderation.R` and `test-kit243-summary-current-outputs.R` (the section-9 grep at the
  base gives 71 including these two);
- KIT-236: the 34 files of its correction run (`kit236-c1-8fd559b2`);
- 10 files are shared (each member runs the other's: `test-e2e-browser-harness.R`, `test-e2e-proteomics-dia-limpa.R`,
  `test-e2e-proteomics-tmt-lfq.R`, `test-enrichment-significance.R`, `test-gui-entry-proteomics.R`,
  `test-kit205-enrichment-background.R`, `test-kit215-not-tested-lists.R`, `test-kit220-da-threshold-consumers.R`,
  `test-module-ci-prot-da.R`, `test-omic-summary-module.R`); KIT-249's new file is in KIT-236's scope only by this
  union, and KIT-236's new file in KIT-249's.
- Script on the head: `docs/qa/workflow-execution/kit/kit-236/render.R`, no arguments, out `kit236-render`, pattern
  `\.(txt|tsv)$`, timeout 90 minutes, exactly as the correction run 45621749264 requested it.

**Run settings.** `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` on both sides (`test-prot-norm-core.R` is not in scope, so no
`MULTISCHOLAR_KIT105_W3_DUMP`); no extra packages (neither member's runs used any); reference BLAS; 8 shards per side;
per-file kill limits of 120 minutes for the three GUI-entry files and 90 for `test-workflow-saved-da.R` (not
thresholds).

**Ceiling.** 1,560 job-minutes: KIT-249's merge-check (R2) ceiling 1,020 (12 × 70 + 2 × 90) and KIT-236's R2 ceiling
540, each from its ticket's budget.

**Base values.** Per file in `merge-train-20261126e-predictions.tsv`. Train 11 moves no file in this scope (its
declared moves are `test-prot-norm-core.R` and the KIT-257 file, neither here), so base counts are staging's: train 11's
declared table (the checkpoint 45644467847 at 828dd105), else train 10's, 9's or 8's run tables. The checkpoint
installed `sva`, `limpa` and `mixOmics` and train 8 `sva`; this run installs none, so skip counts in files that need
them may read differently; the verdict is base against head in this run. Nine files are on no record
(`test-prot-norm-optimization-direct-shared.R`, `test-protein-da-golden-master.R`, `test-prot-ruv-sweep-imputation.R`,
`test-lipid-12c-ruv-s4-characterization.R`, `test-metab-02y-get-neg-ctrl-metab-anova-characterization.R`,
`test-e2e-manifest.R`, `test-module-ci-prot-import.R`, `test-module-ci-sentinels.R`,
`test-prot-summary-report-templates.R`): predicted `same`, counts read from the run.

### Predictions, KIT-249

- `test-kit249-prot-da-qvalue-fallback.R`: only on head, 144/0.
- `test-workflow-saved-da.R`: `message-changed`, base 776/175, head 755/196. Newly failing on head exactly R2's 21
  info lines (`kit-249/r2-45601895895/saved-da-new.txt`: 18 `entities`, 3 `calls`); the messages changed are exactly
  I01's 87 leaves; KIT-201's four `record` leaves, `limpa_one_contrast`, `states`, `fits` and `fdr` do not move; the
  digest values in the changed messages are read from the run. Nothing else in the file moves.
- `test-kit201-treat-moderation.R` 89/0 and `test-kit243-summary-current-outputs.R` 246/0: `same`.
- Every other file of its 79 `same`.

### Predictions, KIT-236

- `test-kit236-report-sentences.R`: only on head, 289/0.
- `test-kit227-report-ruv-sentences.R`: head 162/0; base 166/0 (train 4's run 45586240983; the file and the
  lipid/metab templates it reads are unchanged since); `count-changed` by KIT-236's declared edit; base count read
  from the run.
- `test-kit235-lipid-metab-report-parser.R`: head 54/0; base 50/0 (train 4); `count-changed`; base count read from
  the run.
- The other 31 `same`, with the existing failure lists identical: `test-gui-entry-proteomics.R` 168/25,
  `test-omics-summary-core.R` 676/25, `test-omic-summary-module.R` 1675/28, `test-kit220-da-threshold-consumers.R`
  79/0 (with the trx flag; KIT-236's own run, without it, read 71/2).
- Render job: exit 0, `checks.tsv` 140 rows, 140 "as predicted".

### Combined predictions, shared files

The ten shared files are `same` under both members at the base counts in the table (for example
`test-gui-entry-proteomics.R` 168/25, `test-kit205-enrichment-background.R` 228/0, `test-module-ci-prot-da.R`
250/0, `test-omic-summary-module.R` 1675/28). KIT-249 changes a proteomics comparison's q-values only when Storey's
q-value fails; no GUI-entry lane or module-CI fixture in these files has such a comparison (KIT-249 R2: the proteomics
GUI-entry oracles identical; KIT-236 closure: no lane has a fallback). KIT-236 changes report templates only; these
files read no template sentence it changes (its own run: all `same`). `test-workflow-saved-da.R` is KIT-249's alone
(KIT-236 does not run it).

### Predictions, R CMD check

- Findings identical on both sides, nothing only on head; the portable-file-names NOTE the same count on both sides.

### Decision and attempts

- Every prediction holds: TRAIN GREEN, landing after train 11 lands with its code unchanged.
- Any other move is a result: the member is named with its files and counts; it is not repaired here.
- One run; one re-send only for an instrument defect with no result judged. An HTTP 500 at dispatch with no run
  created: started once by hand with the README's command, recorded.

## Result

Run 45667562462 (request `merge-train-20261126e-55413b59`, run-tests commit 955c0722, dispatched normally), head
55413b59 against train 11's tip a50f9302, sent once; conclusion success; 565.9 job-minutes of the 1,560 ceiling (tests
518.1, check 26.1, render job 20.8, plan and results 1.0). Tables in `merge-train-20261126e/` (`compare.tsv`,
`tests.tsv`, `tests-base.tsv`, `check-compare.txt`, `check-summary.txt`, `meta.txt`, `render-checks.tsv`).

Observed, 103 files: 98 `same`, 5 not `same`:

| file | base p/f/e/s | head p/f/e/s | verdict | prediction |
|---|---|---|---|---|
| test-kit249-prot-da-qvalue-fallback.R | absent | 144/0/0/0 | only-on-head | KIT-249, held |
| test-workflow-saved-da.R | 776/175/0/0 | 755/196/0/0 | message-changed | KIT-249: counts held; newly failing held; changed messages 81, NOT the predicted 87 (below) |
| test-kit236-report-sentences.R | absent | 289/0/0/0 | only-on-head | KIT-236, held |
| test-kit227-report-ruv-sentences.R | 166/0/0/0 | 162/0/0/0 | count-changed | KIT-236, held |
| test-kit235-lipid-metab-report-parser.R | 50/0/0/0 | 54/0/0/0 | count-changed | KIT-236, held |

Every `same` file has the base counts in `merge-train-20261126e-predictions.tsv` where they were on record (0
differences, script over `compare.tsv`), including the ten shared files (`test-gui-entry-proteomics.R` 168/25,
`test-kit220-da-threshold-consumers.R` 79/0, `test-omic-summary-module.R` 1675/28, ...), `test-kit201-treat-moderation.R`
89/0 and `test-kit243-summary-current-outputs.R` 246/0. Read from the run for the nine files not on record: see
`compare.tsv` (all `same`).

`test-workflow-saved-da.R`, parsed from `failures-base.txt` and `failures.txt` (each failure keyed by test and leaf; no
key repeats on either side):
- 21 newly failing on head, 0 only on base: the set equals R2's `saved-da-new.txt` exactly (held).
- 175 failing on both sides: 94 identical, 81 with changed messages. Predicted 87 (I01's table). The 81 are
  `products` 18, `dependencies` 18, `head_manifest` 12, `display` 11, `manifest` 11, `files` 11.
- The 6 leaves of I01's 87 that did NOT change: the `files` leaf of the three qvalue cases (dia_changed_arguments,
  dia_qvalue_only, dia_refused), in W1 (a) and in W2. Each fails on both sides with an identical message (for example
  dia_qvalue_only: actual abadbba0…, expected f5079a9e…).
- No `record`, `states`, `fits`, `fdr` or `limpa_one_contrast` leaf moved (held).
- Inferred, not shown: I01 attributed the qvalue cases' `files` move to earlier clicks' files left in the shared
  publication folder. On this base KIT-243 (staging since train 8) makes each DA remove its own earlier files, and
  train 8 recorded these three cases' `files` messages as KIT-243's; so KIT-249's change to the fallback clicks' files
  may no longer reach these three digests. The register's prediction was by reading on a base without KIT-243 in that
  form; this is recorded as a deviation, not as held.

Render job (KIT-236): exit 0; `checks.tsv` 140 rows, 140 "as predicted"; console "0 of 140 checks differ from the
prediction" (held).

R CMD check: 0 ERRORs, 3 WARNINGs, 3 NOTEs on both sides; all 6 findings identical, none only on head or only on base.

Verdicts:
- KIT-236: holds exactly (new file 289/0; KIT-227 166 -> 162, KIT-235 50 -> 54; the other 31 `same`; render 140/140).
- KIT-249: holds except one count: `test-workflow-saved-da.R` changed messages 81 against the declared 87 (6 qvalue-case
  `files` leaves unchanged). No move outside the declared set; the deviation is fewer moves, not more.
- Check: identical.

Train verdict: NOT GREEN as declared, on that one KIT-249 count; every other prediction held. The coordinator decides
(the stop rule named "any other move"; this is a declared move that did not happen). The train is not rebuilt here. It
lands, if at all, only after train 11 lands with its code unchanged.

Not checked: train 11's own result; any file outside the 103; the digest values (read from the run, captured at the
re-record under VAL-274); the full suite (each member's A4 at the next checkpoint).

## Merge record (coordinator)

Landed on `v0.9-staging` on 2026-11-26 by a merge of `merge-train-20261126e` (tip c19f02bd) into
staging 12701142. The stacking condition held: train 11 landed with its code unchanged (0a774e0b,
KIT-257 = VAL-281), and staging differs from the commit this train was built on (a50f9302) only
under `docs/`. The only conflict at landing was the ledger: staging's rows kept through VAL-294,
then this train's rows. Ledger rows: KIT-249 = VAL-295; the coordinator's correction of VAL-237's
evidence sentence = VAL-296 (the KIT-249 audit's text); KIT-236 = VAL-297.

The deviation and its ruling. The declaration predicted 87 changed messages in
`tests/testthat/test-workflow-saved-da.R`; the run showed 81. The six that did not move are the
`files` leaf of the three TREAT-0.5 cases (dia_changed_arguments, dia_qvalue_only, dia_refused) in
W1 (a) and W2, identical on both sides. KIT-249's auditor settled it from the failure records and
the code ("Train 12 result", AUDIT-KIT-249-impl-20261125): on a base with KIT-243 each analysis
removes its own earlier files before it publishes, so those listings no longer carry the fallback
clicks' output; the prediction of 87 was the auditor's error (true only on a base without
KIT-243); the fix is shown acting in every fallback click (a page saved where q-values were empty
before). The declared stop rule concerned moves outside the declared set; none occurred. Landed on
that ruling; KIT-249's ledger row states both counts.

Evidence: run 45667562462 (base a50f9302, head 55413b59; 103 test files, R CMD check and KIT-236's
render job): 98 files the same; the two new test files only on head and passing (144/0 and 289/0);
two declared count changes; `test-workflow-saved-da.R` 776/175 to 755/196 with exactly the 21
declared newly failing expectations; the render job 140 of 140 checks as predicted; the check
identical on both sides (0 errors, 3 warnings, 3 notes). 565.9 of 1,560 job-minutes.

Carried by this landing: the KIT-112 captures of `test-workflow-saved-da.R` wait for the
end-of-project re-record (VAL-274) and are not to be read as correct output until then; KIT-236's
carried items for the second of KIT-236 and KIT-242 (the per-assay RUV factor counts, the
normalisation and internal-standard sentences, the caption's method words, the volcano sentence
at cutoffs above 0). Open with the owner: whether the notice names the comparisons that fell back
to Benjamini-Hochberg and a sentence that BH is more conservative than Storey's q.

Not shown by this record: the full suite on the landed tree (next checkpoint).
