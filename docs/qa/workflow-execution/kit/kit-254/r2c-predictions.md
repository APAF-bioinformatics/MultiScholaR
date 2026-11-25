# KIT-254 correction pass and follow-up R2c: predictions, written before the request

Base 44e3e910 (kit-254: the fix, staging c4505c7f merged, the correction pass's tests, no statement code); head
904a5514 (the statements, `R/utils_workflow_analysis_record.R` only). 58 files: R2's scope regenerated at c4505c7f (48,
with `test-kit224-ruv-feature-levels.R`), `test-kit214-sample-columns.R`, and the ten replay files (the nine consumers
and `test-workflow-analysis-record.R`); longest first, 4 shards per side, R CMD check on both.

- `test-kit254-omics-correlation-skip.R`: `fixed`: base 20 failures (E8 4 per omic, E15 and E16 2 each per omic, E17
  1 per omic, E18 2), head 0; every other row passes on both (G9 to G13 included).
- `test-omics-norm-state-guards.R`: `same` (651/0 both sides: the T1/T8 edit is on both).
- Every other file `same`, including the seven with pre-existing failures (I05) at staging's counts and the replay
  consumers: none of their records has a lipid/metab selection followed by the skip step or a selection whose input is
  not the step before it, so no statement and no warning is added (by reading).
- R CMD check: findings identical, none only on head.
