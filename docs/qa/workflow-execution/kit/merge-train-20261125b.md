# Merge train 2 of 2026-11-25: KIT-200

Coordinator's instruction after train 1 landed (staging 185ebf75): same procedure as train 1
(`merge-train-20261125.md`).

## Members

- KIT-200 (`origin/kit-200-prot-da-proteins-without-values` at 269a7f7a): implementation audit
  and closure check, AUDIT-KIT-200-impl-20261125, "ready to merge by the next train"; I01 to I07
  closed. Three non-blocking text corrections C-a to C-c in the register are not applied by the
  train (records only; left for the coordinator).
- Not members: KIT-227 (merges after KIT-220, which is not ready); KIT-215 (implementation audit
  running). No other branch is recorded in `tickets/kit/audits/DISPATCH-STATE.md` as audited and
  ready to merge with no open repair.

## Declaration (written before the run is sent)

- Purpose: KIT-200's check on the merged tree (closure check, "For the next train").
- Base commit: `185ebf75cc179fdf88df397c38473cdf9592bb6f` (`origin/v0.9-staging`, train 1's head).
- Head commit (bound by the run): `6f206314` (full hash in the run request), the merge of KIT-200.
  This file and the records commit after the run are under `docs/qa/` only.
- Conflicts met: `docs/methods-decision-ledger.md` only (staging's rows to VAL-200 kept in order,
  KIT-200's row after them written `VAL-<train2-1>` for the coordinator to number when landing
  the train, since owner-decision rows VAL-201 onward are being added to staging). NEWS merged
  without conflict, under the one unreleased heading; `DESCRIPTION` unchanged. Code changed
  relative to the base: `R/func_prot_da_model_stats.R`, `man/runTestsContrasts.Rd`, the new test
  file.
- Scope: KIT-200's R2b list (47 files, `kit-200/r2b/compare.tsv`). The ticket's section-9 grep
  on the train tree adds only three helper files (`helper-*`, dropped by the ticket's rule), so
  nothing is added. R CMD check on both sides. Settings as R2b (3 shards per side, no extra
  packages, no env), except the per-file limit of `test-gui-entry-proteomics.R` at 120 minutes
  (as train 1; a kill limit, not a threshold).
- Ceiling: 540 job-minutes (R2b's 36.0 core-hours).
- Decision: every file `same`, apart from `test-kit200-prot-da-no-value-proteins.R` only on head
  and passing (538 in R2b); `test-gui-entry-proteomics.R` `same`; R CMD check findings identical
  on both sides or fewer on head. Judged against train-1 staging (closure check). If it holds:
  TRAIN GREEN. If not: the files and counts are reported; the train is not rebuilt here.
- Attempts: one run; one re-send only for an instrument defect with no result judged.

## Result

(filled after the run)
