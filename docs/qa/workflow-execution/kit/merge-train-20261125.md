# Merge train 2026-11-25: KIT-202 then KIT-184

Coordinator's ruling of 2026-11-25 (`tickets/kit/audits/DISPATCH-STATE.md`, "merge train"):
the two audited tickets merge as one train, checked by one run on the train's head against
staging as base. KIT-200 is not in this train (it joins a later one).

## Declaration (written before the run is sent)

- Purpose: the check on the merged tree that each ticket's closure requires (KIT-202: R2b on
  the merge, AUDIT-KIT-202-impl-20261125 closure check; KIT-184: the merge-condition check,
  AUDIT-KIT-184-impl-20261125 Re-check 1), done once for both on the train.
- Base commit: `cf8af2453a3766c319f0aec42c1bde8b4cc68ed2` (`origin/v0.9-staging` when the train
  was built).
- Head commit (bound by the run): `485ec43d` (full hash in the run request), the second merge.
  Train commits:
  - `01929e1a` Merge KIT-202 (branch at 6370d5b3); ledger row VAL-197.
  - `485ec43d` Merge KIT-184 (branch at 7920a818); ledger rows VAL-198, VAL-199, VAL-200.
  This file and the records commit after the run are under `docs/qa/` only.
- Conflicts met: one per merge, `docs/methods-decision-ledger.md` only (staging's rows kept in
  order, the ticket's rows after them, numbered with the next free numbers). NEWS merged without
  conflict; both tickets' entries are under the one unreleased heading; `DESCRIPTION` unchanged.
- Scope: the union of KIT-202's R2 list (43 files, request `kit202-r2-6e02de15`) and KIT-184's
  R2 list as last run on a merge (29 files, request `k184-r3b-d7928ac9`, which includes
  `test-gui-entry-proteomics.R` with the GUI-entry lanes): 67 files (5 shared), base and head,
  R CMD check on both sides. `sva` and `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` as KIT-184's runs
  (both sides alike). Per-file limits 120 minutes for the four slow files as before.
  No scripts (the request's `timeout_minutes` field belongs to scripts only).
- Ceiling: 900 job-minutes.
- Decision: every file `same` or better, apart from each ticket's declared differences:
  - KIT-202: `test-kit202-intensity-filter-samples.R` only on head and passing;
    `test-metab-qc-shared.R` count-changed with no failure; `test-workflow-omics-declarations.R`
    and `test-workflow-omics-binding-steps.R` message-changed (VAL-181, VAL-189).
  - KIT-184: `test-integration-removed-factors.R` only on head and passing;
    `test-gui-entry-proteomics.R` `same` (lane `prot_dia`'s 4 failures on both sides;
    `prot_dia_limpa` passing on its re-frozen reference).
  - R CMD check findings identical on both sides, or fewer on head.
  If all hold: TRAIN GREEN, the coordinator may fast-forward staging to the train's head
  (provided staging is still at the base). If not: the failing ticket is named with its files
  and counts; the train is not rebuilt here.
- Attempts: one run. A failure in the instrument's mechanics with no result judged may be fixed
  and sent once more (recorded here); a second failure stops the train.

## Result

Run 45503677179 (request `merge-train-20261125-485ec43d`, run-tests commit 664d0579), head
485ec43d against base cf8af245, sent once; conclusion success; 538.7 job-minutes of the 900
ceiling (tests 511.6, check 26.2, plan and results 0.8). Tables in `merge-train-20261125/`
(`compare.tsv`, `tests.tsv`, `tests-base.tsv`, `check-compare.txt`, `check-summary.txt`,
`meta.txt`). Staging was still at cf8af245 when the result was read.

Observed, 67 files: 54 `same`, 13 not `same`:

| file | base p/f | head p/f | verdict | ticket |
|---|---|---|---|---|
| test-kit202-intensity-filter-samples.R | absent | 163/0 | only-on-head | KIT-202, declared |
| test-metab-qc-shared.R | 19/0 | 20/0 | count-changed | KIT-202, declared |
| test-workflow-omics-declarations.R | 46/22 | 45/23 | message-changed | KIT-202, declared (VAL-181) |
| test-workflow-omics-binding-steps.R | 58/1 | 58/1 | message-changed | KIT-202, declared (VAL-189) |
| test-integration-removed-factors.R | absent | 526/0 (3 skipped) | only-on-head | KIT-184, declared |
| test-workflow-integration-input.R | 181/0 | 241/0 | count-changed | KIT-184 |
| test-integration-input-lipid-metab.R | 132/0 | 208/0 | count-changed | KIT-184 |
| test-integration-input-transcriptomics.R | 57/0 | 78/0 | count-changed | KIT-184 |
| test-integration-enrichment-tests.R | 180/0 | 181/0 | count-changed | KIT-184 |
| test-integration-gui.R | 395/0 | 397/0 | count-changed | KIT-184 |
| test-integration-enrichment-mofa.R | 0/0 (2 skipped) | 0/0 (3 skipped) | count-changed | KIT-184 |
| test-omic-summary-module.R | 1675/28 | 1675/28 | message-changed | KIT-184 |
| test-omics-summary-core.R | 676/25 | 676/25 | message-changed | KIT-184 |

`test-gui-entry-proteomics.R` is `same` (189/4 on both sides; the 4 are lane `prot_dia`);
`test-gui-entry-lipidomics.R` and `test-gui-entry-metabolomics.R` `same`. R CMD check: 0 ERRORs,
3 WARNINGs, 3 NOTEs on both sides, all 6 findings identical, none only on head.

Each of the four KIT-202 rows above is byte-identical to its row in KIT-202's R3 table
(`kit-202/r3/compare.tsv`), and each of the nine KIT-184 rows is byte-identical to its row in
KIT-184's audited R3b table (`kit-184/r3/b/compare.tsv`).

Verdicts:

- KIT-202: meets its declared differences exactly; every other file `same`.
- KIT-184: its two declared items met (new file only on head and passing; GUI-entry file `same`).
  Eight further rows are not in the coordinator's short list: six count-changed files where
  head passes more tests and fails none (this ticket's added tests), and the two summary files,
  message-changed on expectations that fail on both sides. In those two the transcriptomics
  captures now list the integration file the ticket writes ("transcriptomics": "<present>").
  This is the ticket's intended change against references that already fail on staging,
  recorded in `kit-184/outcome.md` and accepted at Re-check 1 on R3b. Inferred: same or better
  under the ticket's audited declarations; the coordinator confirms.
- Check: identical on both sides.

Not checked: the trees after the run's head (docs only: this file, the tables, the ledger's
status cells); any file outside the 67.
