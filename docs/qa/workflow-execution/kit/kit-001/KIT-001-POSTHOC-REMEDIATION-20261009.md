# KIT-001 post-hoc audit: remediation record

Register: [KIT-001-POSTHOC-AUDIT-20261009.md](KIT-001-POSTHOC-AUDIT-20261009.md) (not rewritten).

| Finding | State | Closure evidence |
|---|---|---|
| PA-08 | closed 2026-10-09 | User decision: the ceiling removal covers every proteomics intensity-filter site (five); lipidomics and metabolomics keep theirs until KIT-031/032. Plan decision row 2026-10-09. |
| PA-01 | closed 2026-10-09 | successor script `attempt4-post-hoc-correction-v2.R` records the full correction (fingerprints, verifier rerun, decision) as it ran; the in-run fingerprint requirement is now a KIT-005 acceptance item. |
| PA-02 | closed 2026-10-09 | provenance note in `outcome.md`: the correction ran with the content later committed as `dfd6aacd`; the member's recorded `3512f84d` is explained. |
| PA-03 | closed 2026-10-09 | claim scope stated in `outcome.md` (stop between steps); mid-step fault is now a KIT-005 acceptance item. |
| PA-04 | closed 2026-10-09 | install sequence (`sLJ4aC` failed, `uvS55B` superseded, `y8rDOi` bound) recorded in `outcome.md`. |
| PA-05 | closed 2026-10-09 | `kitJourneySteps()` import step state is `initial`; test pins it; harness hash rebound. |
| PA-10 | closed 2026-10-09 | ticket and outcome wording corrected (56 assertions; harness revision of the failure member stated). |
| PA-06, PA-07, PA-09 | informational | no action; PA-06 is first exercised by KIT-005's control journey. |
