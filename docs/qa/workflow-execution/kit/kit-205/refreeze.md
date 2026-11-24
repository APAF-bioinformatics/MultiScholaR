# KIT-205 re-freeze (VAL-174)

The owner allowed the re-freeze as an exception by file (methods ledger VAL-174, staging 9f3b7c49; audit
`AUDIT-KIT-205-impl-20261125.md`, re-check 1). Branch `kit-205-refreeze` from `v0.9-staging` 9f3b7c49.

**The edit (60ebb21f), `tests/testthat/test-omic-enrich-module.R` only.**
- `omicEnrichW1FrozenUi`, lines 138 to 206, replaced by the first 69 lines of `dumps/r2b/dump-ui-head-1.R`; the end
  marker at 207 stays.
- The five `omicEnrichW113Frozen` cases `gprofiler`, `kernel_gprofiler`, `restore`, `organism_filter/no_mapping`,
  `organism_filter/mapped` replaced by the same-named blocks of `dumps/r2b/dump-w113-head-1.R` (each one line longer).
- Diff (19 insertions, 14 deletions): the UI block's HTML digest and the three texts "measured" to "tested"; per W113
  case the provenance `rest` digest and one new record line `parameters$enrichment_background` = "fixed ; fixed ;
  tested". Nothing else in the file changed; the two clusterProfiler W113 cases and every `omicEnrichW1Frozen` capture
  are untouched. Control: KIT-205 R2b, run 45404531700 (`dumps/r2b/comparison.md`).

**Closure run 45465468140** (declared before sending: `test-omic-enrich-module.R`, base 9f3b7c49 against head 60ebb21f,
one shard per side, ceiling 120 job-minutes; 36.3 job-minutes, 2.42 core-hours). Result (`refreeze/`): the file goes
from 1042 passed / 28 failed to 1051 / 19.
- Test (a) `:209`: `fixed` (failing on base, passing on head).
- R3 P-28 `:4309` (the comparison at `:4321`): `fixed`.
- The W113 test `:4201`: `message-changed`; on head it fails only for `clusterprofiler` and `kernel_clusterprofiler`
  (KIT-173's plan section 11 row); the five re-frozen cases pass.
- Everything else in the file is unchanged, including the `(b) clusterprofiler` and `(b) kernel_clusterprofiler`
  captures (failing on both sides). The `gprofiler_plot` flake did not show.

As expected; no other difference.
