# KIT-254 follow-up R2 (OD-1 (a), VAL-247): predictions, written before the request

Base ba383de2 (the follow-up's tests-only commit on kit-254 merged with staging c4505c7f); head = the statement
commit. Scope: 57 files at the merged tree (R2's 46 existing files, KIT-224's `test-kit224-ruv-feature-levels.R`
which the grep now finds, the nine replay consumers of ticket section 10, and the new file), longest first, 4 shards per
side, R CMD check on both.

- `test-kit254-omics-correlation-skip.R`: `fixed` (base 8 failures, the E8-statement rows; head 0).
- Every other file `same`, including `test-workflow-analysis-record.R` and the replay consumers: none of their records
  has a lipid/metab `<prefix>_norm_complete` step directly after `<prefix>_correlation_filtered` (by reading: KIT-212
  G5 makes such a record but does not replay it), so no statement and no warning is added to any of their replays.
- Check: findings identical on both sides, none only on head (`@noRd` functions; `man/`, `NAMESPACE` unchanged).
