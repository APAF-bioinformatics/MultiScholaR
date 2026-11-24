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

(filled after the run)
