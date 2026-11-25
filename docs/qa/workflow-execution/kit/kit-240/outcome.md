# KIT-240 outcome: the integration tabs on by default

Ticket revision 1.6. Branch `kit-240` from `v0.9-staging` 3eeaf100. Every run on `run-tests` (R 4.6.1, reference
BLAS, ubuntu-22.04); run ids, commits and job minutes in `runs.tsv`. Evidence copies under `r1/`, `r2/`, `r3/`, `r4/`, `r5/`, `merged/` (I and the base-against-head run), `m/` (M) and `b/` (B).

| Commit | Content |
|---|---|
| 0470164c | tests only: `test-kit240-integration-default.R` (T1, T3, T4, T6, T7, T8, T9); T2 in `test-integration-gui.R` |
| 8cdd1f6e | fix: package default TRUE; `entry.R` and `lifecycle-test.js` (RUV panel only); `missing` in the probe and the Model tab's sentences; marks reworded (values unchanged); `elbo_last`, `elbo_rank`, the Results and record sentences; moved tests; README, NEWS, guide, ledger `VAL-269`; `default-identity.R`, `wording.md` |
| ba6b44d8 | tests only (r1.3): T10, the Enrichment tab's mark |
| fac10b51 | fix (r1.3): `integrationEnrichmentLabel` in `R/mod_integration_enrichment.R`; `mark_provisional` "default not yet confirmed"; the identity script's allowed insertion |
| b89be231 | tests only (r1.4): T7 (I01), T8 (I09), T12 (I03, I05, I09, I10), `test-integration-gui.R` (I02) |
| 37071446 | wording and records (r1.4): I01 display, I02, I03, I04, I05, I06 NEWS, I07, I09, I10 |
| 651077e7 | records: R4 |
| 925e532e | tests only (r1.5): the VAL-252 clause in T7 and T8 |
| c74e3b36 | the VAL-252 clause on the Results tab and in the Methods paragraph (draft) |
| ccc0a401 | records: R5 |
| bded5eb2 | merge of staging 63143aa3 (KIT-234, KIT-239): ledger, NEWS, the Results list resolved as AUDIT-KIT-240-impl section 6 |
| 83d1b222 | records: I, M and B at bded5eb2 |
| cf91f07a | C1 (closure check): the clause's two strings corrected (draft) and `wording.md` |
| (this commit) | C2 records; then the merge of staging 42846c0e |

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

R5 (45603743811; base 925e532e, head c74e3b36; evidence under `r5/`): base 2 new failures, 0 errors, both in T7 (the
tab clause and the Methods clause), as declared; head `test-kit240-integration-default.R` 92 passed, 0 failed (`fixed`,
T9 skipped); `test-integration-results.R`, `-record.R`, `-outputs.R` `same`, 0 failures. Identity script against
3eeaf100: every check TRUE, counts 1, 2, 21 and 8 per mode (no new allowed substitution: the clause is server-rendered).
78.0 job-minutes, 5.20 core-hours (declared ceiling 6.00).

## After merging staging (r1.6): runs I, M and B at bded5eb2

Why these runs did not wait for KIT-228 (ticket r1.6): the ticket names KIT-228 for its content and merge order only;
it changes the omic Summary Save, nothing M, B or I observe. Staging 63143aa3 (KIT-234, KIT-239) was merged into
`kit-240` as bded5eb2, with the conflicts resolved as AUDIT-KIT-240-impl section 6 directs (ledger: staging's rows
VAL-226 to VAL-259 kept, KIT-240's row after them; NEWS: both entries; `R/utils_integration_results.R`: `limitations`
then `elbo_rank`; KIT-239's `help_scaling` exactly as on staging, VAL-251). KIT-228 merges second as a partner and
re-runs what it touches, `test-kit240-integration-default.R` and `kit-240/default-identity.R` included.

- **Base against head and I (45608060716, base 63143aa3, head bded5eb2; `merged/`).** As declared: 0 `new-failure`,
  0 `message-changed`; the new file `only-on-head` 92 passed, 0 failed, 1 skipped (T9); `test-integration-gui.R`
  401 -> 405 and `test-integration-model.R` 263 -> 265 (`count-changed`, the moved and added assertions); 24 files
  `same`, `test-workflow-checkpoint-harness.R`'s one error identical on both sides. R CMD check 0 ERROR, 3 WARNING,
  3 NOTE on both, 0 only on head. Identity against the merge base: every check TRUE; per mode settings note 1, Results
  note 2, mark 21, Enrichment insertion 8 on head and 0 on base. **A2 met.**
- **M (45608049489, `run-mofa-smoke-k240`; `m/`).** MOFA2 1.22.0, reticulate 1.46.0, mofapy2 0.7.3. Six files, 0
  failures; one skip, T4 (by design: all software present). The reported model's factors, weights and variance
  explained equal `mofa_reference/` within 1e-5 (`test-integration-mofa-model.R`, `-model-tab-mofa.R`); the workbook,
  record and Results paths pass with real MOFA2. T9 (`m/t9-elbo-lines.txt`): the five seeds' final ELBOs are
  -45695.190153, -45695.191540, -45695.190312, -45695.191179, -45695.191303; the declared seed ranks 1 of 5; range
  0.0014 (relative 3.0e-8): the first real-corpus evidence for I01 and VAL-252's clause (a spread far below a few
  rounds' change). The record's Methods prints "0.00 below the highest (range over the models less than 0.01", then
  the clause. **A3 met.**
- **B (45608074480, stage full, both platforms; `b/`).** All nine lifecycle passes CHECK PASSED; the tile is live in
  every pass; pass 2's R line reports the variable and `ruv_design_panel` true with no `integration_gui` key; passes 1
  and 3 have no test-features line. The version table check passed (complete, no hard-closure package missing);
  MOFA2 1.22.0 with its mofapy2 0.7.3 pins and `mofapy2 0.7.3` printed; uwot and Rtsne named among the bundled
  packages; reticulate's and clusterProfiler's rows are not printed in the log (covered only by the completeness
  check). 0 artifacts uploaded. The three jobs are red on their GUI-entry lane steps; the lane tables are identical,
  per platform, to run 45386641894 (branch kit-167, 3c4c5c90, before KIT-240): stored laptop oracles of other
  tickets, not KIT-240 (`b/lanes-*.tsv`). **A4 met** (lifecycle); the red lane steps are recorded, not this ticket's.
- Job minutes: 205.3 (13.69 core-hours), 23.0 (1.53), 245.4 (16.36).

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
- **Evidence.** R1 45579469695, R2 45581871256, R3 45587886914, R4 45596337362, R5 45603743811 (run-tests); at the
  merge with staging 63143aa3 (bded5eb2): I with the base-against-head run 45608060716, M 45608049489 (MOFA lane, T9's
  final ELBOs, `mofa_reference/` within 1e-5), B 45608074480 (bundle lane, three images, nine lifecycle passes);
  `runs.tsv` has commits, minutes and results. The train's run repeats I and the run-tests files; the checkpoint's
  full-suite run id is owed (A7).

## Closure check (2026-11-26): C1 to C4

- **C1** (major, repaired in cf91f07a): the VAL-252 clause said "a difference of this size" / "Differences of this size
  were observed", false for runs whose seeds differ far more (KIT-151c's uncapped arm, `feature_selection = "none"`:
  258.84 below the highest, against a largest observed one-solution difference of 0.676), and it called the KIT-062
  pair "one solution trained twice" without saying it was single precision (the shipped setting is double, where that
  diagnosis found the runs identical). Both strings now carry the audit's texts; drafts for the owner. Checked by
  reading: T7 reads both keys from the wording tables; T8's "where training stopped" is in both texts and they print
  only with a rank; no other test and not the identity script pins them (they are server-rendered). The train's run is
  the check.
- **C2** (records): this commit (evidence folders, commit table, the change report's evidence, the ledger row,
  `wording.md`'s revision, B's step conclusions `b/steps-45608074480.tsv` and version-table lines
  `b/version-table-summary.txt`, the KIT-228 note in its ticket).
- **B's steps** (`b/steps-45608074480.tsv`): in each of the three jobs the only failed step is the GUI-entry lanes
  step; every later step ran (they depend on the install and test tools, not on the lanes step) and passed or was
  skipped by its own condition (parity, D5, the 24.04 e2e lanes by matrix, a post step). Version table
  (`b/version-table-summary.txt`, all three jobs): complete, no hard-closure package missing; "in the lockfile, not
  bundled" names neither reticulate nor clusterProfiler, so both are bundled on each image; uwot and Rtsne are named as
  bundled.
- **C3** (note, the owner's polish): at rank 1 the tab and the Methods print "ranks 1 of 5: 0.00 below the highest",
  true but redundant beside the rank.
- **C4** (note): `methods_provisional` still says "provisional defaults" while the status word is "default not yet
  confirmed"; no Model or Results setting can reach it today.

## Release checklist (owed by this ticket; no checklist file exists yet)

1. **KIT-228 merged before the release commit** (VAL-201: a Save that writes no new integration input removes the
   earlier one), since the Integration tabs now reach every user. KIT-228 merges after KIT-240 and re-runs
   `test-kit240-integration-default.R` and `kit-240/default-identity.R`.
2. **The systemd observation** (R1-05): integration observed once on the systemd worker of a Linux desktop install of
   the release candidate (staging with KIT-240 merged), started without `MULTISCHOLAR_TEST_FEATURES` (also not from
   `~/.profile`; no "(test features on)" title, no `MSR-TEST-FEATURES` line), with the R log's `MSR-LAUNCHER` line saying
   `systemd` and a worker unit's log for the Train in the run folder's systemd launcher folder. The owner's session
   (KIT-170 6.8), not agent compute.
3. **The eight enrichment defaults** (R1-04): confirmed by the owner before the release commit, or the release notes
   state that the Enrichment tab's defaults (Pearson correlation, at least 5 samples, STRING species 9606, KEGG hsa,
   Reactome Homo sapiens, pathway sizes 5 to 500, seed 20261121) await confirmation.

## Not checked

- KIT-228's own merge check (it merges second and re-runs this ticket's test file and identity script).
- reticulate's and clusterProfiler's version rows in B (only the completeness check covers them).
- The GUI-entry lanes' stored-oracle failures (identical to before KIT-240; other tickets').
- The systemd controller showing the tabs (reading only; release-checklist item).
- A7 (checkpoint full suite), A5 (owner's wording).
