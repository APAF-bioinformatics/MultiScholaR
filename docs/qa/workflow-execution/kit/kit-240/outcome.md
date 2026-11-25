# KIT-240 outcome (branch point; post-rebase runs owed): the integration tabs on by default

Ticket revision 1.4. Branch `kit-240` from `v0.9-staging` 3eeaf100. Every run on `run-tests` (R 4.6.1, reference
BLAS, ubuntu-22.04); run ids, commits and job minutes in `runs.tsv`. Evidence copies under `r1/` and `r2/`.

| Commit | Content |
|---|---|
| 0470164c | tests only: `test-kit240-integration-default.R` (T1, T3, T4, T6, T7, T8, T9); T2 in `test-integration-gui.R` |
| 8cdd1f6e | fix: package default TRUE; `entry.R` and `lifecycle-test.js` (RUV panel only); `missing` in the probe and the Model tab's sentences; marks reworded (values unchanged); `elbo_last`, `elbo_rank`, the Results and record sentences; moved tests; README, NEWS, guide, ledger `VAL-<merge>`; `default-identity.R`, `wording.md` |
| ba6b44d8 | tests only (r1.3): T10, the Enrichment tab's mark |
| fac10b51 | fix (r1.3): `integrationEnrichmentLabel` in `R/mod_integration_enrichment.R`; `mark_provisional` "default not yet confirmed"; the identity script's allowed insertion |
| b89be231 | tests only (r1.4): T7 (I01), T8 (I09), T12 (I03, I05, I09, I10), `test-integration-gui.R` (I02) |
| (next) | wording and records (r1.4): I01 display, I02, I03, I04, I05, I06 NEWS, I07, I09, I10 |

## Observed

- **R1 (A1), 45579469695 at 0470164c.** T1 (option FALSE, disabled tile, the unavailable panel), T3 (the general
  sentence for MOFA2, reticulate and mofapy2), T4 (`missing` NULL; the runner has MOFA2 and reticulate, so the expected
  piece is mofapy2), T6 (a), (b), (d) and T7 (`elbo_last` and `elbo_rank` NULL, no sentences) failed at their base
  values: 36 expectation failures, 0 errors. T6 (c) and T8 passed; T2 (`test-integration-gui.R`, 397) passed; T9
  skipped (no mofapy2 on run-tests). Stop rule 1 did not fire.
- **R2 (A2 at the branch point), 45581871256, 8cdd1f6e against 3eeaf100.** `r2/compare.tsv`: 21 files `same` (among
  them the moved tests in `test-integration-gui.R`, `-record.R`, `-results.R`, `-outputs.R`, each passing on both
  sides with its own assertions), `test-integration-model.R` `count-changed` (263 to 265 passed, 0 failed), the new
  file `only-on-head` (64 passed, 0 failed, T9 skipped). 0 `new-failure`, 0 `message-changed`.
  `test-workflow-checkpoint-harness.R` has 1 error on both sides, identical (stored laptop reference, not this
  ticket's). R CMD check: 0 ERROR, 3 WARNING, 3 NOTE on both; findings only on head: 0 (`r2/check-compare.txt`).
  Stop rule 2 did not fire.
- **Identity script I, in R2** (`r2/default-identity-verdict.txt`, read from the file): all checks TRUE. (i)
  `app_ui(NULL)` head unset = base TRUE byte for byte in both modes; module equal after the printed substitutions
  (settings note 1, Results note 2, mark "(provisional default)" -> "(decided)" 21 per mode, each found in base and
  absent from head); (ii) head FALSE = base unset for both objects and modes, module the unavailable panel; (iii)
  head identical across two processes. This run is at the branch point; the ticket's A2 evidence is the run after the
  rebase (R1-06), still owed.

## Scope addition r1.3 (coordinator's ruling on R1-04)

The Enrichment tab now marks its eight defaults not yet confirmed (correlation, minimum samples, STRING species, KEGG
organism, Reactome organism, smallest and largest pathway, seed) with " (default not yet confirmed)" after the label.
Reason: release users meet these defaults on the tab, and the ticket's "the mark says so" held only in
`integrationSettings()`. Test T10; the identity script allows exactly this insertion (8 per module render in head,
none in base). Run R3 (45587886914; base ba6b44d8 = the second tests-only commit, head fac10b51; evidence under `r3/`):

- T10 failed on base with 11 expectation failures, 0 errors (`r3/failures-base.txt`). The eight labels carried no
  " (default not yet confirmed)" mark: at ba6b44d8 the mark text is still the r1.2 draft, so the test looked for
  "(default awaiting the owner's decision)" after each label and did not find it. The count was 0, not 8. The mark text
  lacked "not yet confirmed" and contained "owner". On head all of T10 passes; the file moved from 65 passed and
  11 failed to 76 passed and 0 failed (`fixed`); T9 skipped on both.
- The four Enrichment-tab files (`test-integration-enrichment.R`, `-enrichment-options.R`, `-enrichment-tests.R`,
  `-removed-factors.R`) are `same`, with 0 failures on either side.
- Identity script against 3eeaf100 (`r3/default-identity-verdict.txt`): all checks TRUE. In each mode the Enrichment
  mark is inserted 8 times on head and never on base. The substitutions are as in R2: settings note 1, Results note 2,
  mark 21.
- 83.5 job-minutes, 5.57 core-hours (declared ceiling 8.00).

## Revision 1.4: the implementation audit's repairs (AUDIT-KIT-240-impl-20261125)

Applied with the register's texts, as a tests-only commit (b89be231) and a wording-and-records commit; checked by R4
(below).

- **I01 (display only).** `integrationElboText` prints a positive distance below 0.005 as "less than 0.01", never
  "0.00", so the rank line cannot read "0.00 below the highest" at rank 2 or lower; rank 1 still reads "0.00". Nothing
  else about the rank changed. How the rank is presented is the owner's: options (a), (b), (c) or as is, listed in
  `wording.md`.
- **I02 (coordinator's ruling: in this ticket's scope, because the tabs now ship on by default and every release user
  reads it).** The Sample Matching tab said "%d features measured in fewer than two samples are left out" and headed
  its column "(fewer than two samples)". At the defaults the alignment drops features seen in fewer than half of the
  kept samples (`integrationMinObserved`, INT-12), as the record's Methods states. The sentence predates KIT-240: the
  text came with KIT-139 (5d1189ee, 2026-11-20), the threshold with KIT-151 (72d1928b, 2026-11-21); at 3eeaf100 it is
  `R/utils_integration_gui.R:150` and `:156`. Now: "%s: %d features measured in fewer than %d of the %d biological
  samples are left out.", the threshold computed over `nrow(preview$aligned$samples)`. The audit's inferred point (the
  preview's sample table counts the kept samples): `test-integration-gui.R` now asserts 5 rows and a threshold of 3 for
  the small project at the defaults, and `test-integration-partitioned-alignment.R:422` asserts the facts equal between
  the partitioned preview and the whole-object alignment; R4 confirmed it (below). At the tab the
  missing-sample rule is always the default (keep missing), so kept samples equal all samples there; the
  drop-incomplete rule is not reachable from the tab.
- **I03** the alignment note names its three rules; **I05** the report's legend and status word read "a default not
  yet confirmed" / "default not yet confirmed" (no "owner"); **I09** the record states the stability limitation as its
  own sentence whenever two or more seeds were trained; **I10** the Methods sentence defines the ELBO and names its
  units. Tests: T7 (I01), T8 (I09), T12 (I03, I05, I09, I10), `test-integration-gui.R` (I02).
- **I04** the testing guide's introduction and section 6 (text only).
- **I06, observed (not changed here):** R1's log of `test-kit240-integration-default.R` shows, on a runner with MOFA2
  1.22.0 and reticulate 1.46.0 and no `RETICULATE_PYTHON`, the software check (`integrationMofaCapability`, T4) making
  reticulate download uv, CPython 3.13.18 (34.0 MiB) and numpy (15.0 MiB), then report mofapy2 missing. The same
  happens at 3eeaf100 in `test-integration-model.R`. From this merge a `run_app()` user with MOFA2 reaches it on opening
  the Model tab's software check; the installers set `RETICULATE_PYTHON` and are unaffected. NEWS now says so. Whether
  the probe should avoid a managed environment is for the owner and coordinator (stop rule 4 keeps it out).
- **I07** NEWS carries no draft marker (drafts are tracked in `wording.md`); the ledger row cites r1.4, the run ids and
  the Enrichment mark.

### R4 (45596337362; base b89be231, head 37071446; evidence under `r4/`)

- Base: 12 new expectation failures, 0 errors: T7 1 (`integrationElboText(0.004)` gave "0.00"), T8 1 (no limit
  sentence), T12 8, `test-integration-gui.R` 2 (the sparse sentence; the column's "two samples"). **Declared 14; the
  difference is the declaration's arithmetic, not the tests':** T12 has 8 assertions (declared as 9), and the I02
  `expect_identical` of the facts cannot fail on base, because R's `sprintf` ignores the extra arguments and returns the
  old sentence (it passed with a warning). Every assertion that can fail first did. Not re-run.
- Head: `test-kit240-integration-default.R` 88 passed, 0 failed (T9 skipped), `test-integration-gui.R` 401 passed,
  0 failed (`fixed` both); `test-integration-partitioned-alignment.R` (253, its facts-equality check `:422`
  included), `-sample-match.R`, `-record.R` and the four Enrichment-tab files `same`, 0 failures. **I02's inferred
  point is confirmed at the defaults:** the preview's sample table holds the 5 kept samples and the threshold is 3, and
  the partitioned preview's facts equal the whole-object alignment's.
- R CMD check: 0 ERROR, 3 WARNING, 3 NOTE on both; 0 only on head.
- Identity script against 3eeaf100: every check TRUE; per mode settings note 1, Results note 2, mark 21, the Enrichment
  insertion 8 on head and 0 on base (no new allowed substitution: r1.4's texts are all server-rendered).
- 157.0 job-minutes, 10.47 core-hours (declared 10.00: exceeded by 0.47, within R2's correction reserve).

## Revision 1.5: the rank's basis (owner's decision VAL-252)

The ELBO rank is kept, with a clause stating its basis built from two observed facts (KIT-151c's closest cell; KIT-062's
one solution 0.676 apart), on the Results tab (`stability_elbo_slack`) and in the Methods paragraph
(`methods_seeds_elbo_slack`), shown only with the rank. Drafts for the owner's wording pass, with evidence, in
`wording.md`. Tests first (T7 expects both, T8 expects neither without a rank), then the code; run R5 (below, in
`runs.tsv`).

## A8: partner edits and release-checklist items

Done (by the coordinator, read by the audit): KIT-170 6.8 step 2 (`KIT-170.md:664`); KIT-233 `:271` superseded
(`KIT-233.md:376`); KIT-140's open item closed (`KIT-140.md:65`); KIT-238's IJ' without the test setting (its
dependency_note).

Owed (not made here):
- KIT-237 5.4, the `test_features` field;
- KIT-169 `:86`, its NEWS clause and A13 `:249`;
- KIT-181 `:68`, `:110`, A10 and its R1 scripts job;
- KIT-168 8.6 and the kit-167 script (base and head on one side of this merge);
- the next checkpoint's request (`kit-240/default-identity.R`, not kit-138's script);
- the release checklist's two items: the systemd observation on the release candidate without the test setting
  (R1-05's requirements), and the eight enrichment defaults confirmed or disclosed (R1-04).

## For the change report

- **What was wrong, who was affected.** The integration interface was built behind an off-by-default switch; release
  users could not reach it. The installer's test setting showed it, but also showed an unfinished RUV design panel and
  kept Linux desktops off their memory-capped worker. Integration defaults were marked "provisional" or attributed to a
  workbook that does not ship. The Sample Matching tab misstated the sparse-feature rule at the defaults (I02, since
  KIT-151). Only the declared seed's ELBO was kept.
- **What changed.** The tabs are on by default (`options(MultiScholaR.integration_gui = FALSE)` hides them); the test
  setting now switches only the RUV panel (it still keeps the application in the plain launcher); the Model tab names
  the missing MOFA2, reticulate or mofapy2 and how to get it; the marks name the decision behind each default and the
  Enrichment tab marks its eight unconfirmed defaults; the Sample Matching tab states the rule applied; each seed's
  final ELBO is kept and the reported model's rank shown with the stability table's limitation.
- **Code.** At 3eeaf100: `R/utils_integration_gui.R:28-30` (the default), `:150`, `:156` (the sparse texts),
  `packaging/electron/entry.R:14-24`, `R/utils_workflow_integration_model.R:277` (the marks), `:1279-1282`,
  `:1568-1592` (the seeds). At the fix commits 8cdd1f6e, fac10b51 and the r1.4 commit: `integrationGuiEnabled`,
  `integrationMofaCapability`, `integrationSeedStability`, `integrationElboRank`, `integrationElboText`,
  `integrationEnrichmentLabel`, `integrationAlignmentFacts`.
- **Results.** No computed value changes; `mofa_reference/` is not re-recorded. The Methods paragraph gains the rank
  and limitation sentences; the record gains `seed_models$elbo_last` and `model$elbo_rank`.
- **Evidence.** R1, R2, R3, R4 (`runs.tsv`); M, B and the identity run after the rebase are owed; the checkpoint's
  full-suite run id is owed (A7).

## Not checked (so far)

- Runs M and B and the identity run after the rebase (they wait for KIT-234, KIT-239, KIT-228).
- The systemd controller showing the tabs (reading only; release-checklist item).
- A7 (checkpoint full suite), A5 (owner's wording).
