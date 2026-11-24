# KIT-201 outcome (route (a); ticket r1.4)

State: tests-only commit, its failing run, the fix and its passing run done; R0, `m1/precision.tsv`, M1, R2 and R3 wait
for the KIT-196 and KIT-200 merges (their base must contain both: M1's VD20 data sets hold features with no value in any
sample, which stop the DA on a tree without them).

## Commits

| Commit | Content |
|---|---|
| ed56c7a6 | tests first: T-A1 to T-A7, the rewritten legacy-contract tests, `m1-sim.R`, `m1/family.tsv` (K = 990); nothing under R/ |
| 2c399d0c | T-A5's protein table without row names (R1 instrument defect); nothing under R/ |
| fa1f9909 | the fix, `man/runTestsContrasts.Rd`, NEWS, two ledger rows |

## Observed

**R1 (run 45328062374, ed56c7a6, tests only, no base).** `r1/tests.tsv`, `r1/failures.txt`.
- T-A1: lipid and metab `length(unique(s2.prior))` = 1 (stated value); s2.prior differs from the trended robust one;
  `df.prior` length 1.
- T-A2: `length(unique(s2.prior))` = 1 (stated, asserted first); p-values differ from `eBayes(trend, robust)`
  (e.g. 0.1727 against 0.1078).
- T-A3: proteomics at 1 and both S4 methods at 0: `length(unique(s2.prior))` = 1.
- T-A4 (guards): pass.
- T-A5: `workflowReplayRuleMarkers$treat_moderation` is NULL (stated value); then an error in the test's own fixture
  (`column_to_rownames`: the protein table had row names) before any DA ran: an instrument defect (VAL-54), repaired in
  2c399d0c with no change to any expected value, re-run in R1-2.
- T-A6: the second assertion fails with the strings equal (stated).
- T-A7: `fdp` 0.385 (stated); the note lacks "10.1%" and the reworded sentence.
- `test-prot-ruv-factors-in-model.R`: the legacy model record lacks `treat_moderation` (the new field, 1 failure); the
  rewritten legacy test fails at `length(unique(s2.prior))` = 1 and unequal p-values for entry, entry_technical and kernel
  (6 failures). `test-omic-norm-ruv-design.R`: the classic Legacy figure 0.385 against 0.389 (1 failure).
- `m1-sim.R family` on the runner wrote a `family.tsv` identical to the committed one (diff empty).

**R1-2 (run 45329779064; base 2c399d0c, head fa1f9909; one job per side).** `r1-2/`.
- Base (tests-only, T-A5 repaired): T-A5 now runs through and fails at the marker NULL (stated), then at (b)'s record
  and replayed model lacking the marker, "(a) and (b) differ" (base computes one prior under both rules) and (c)'s model,
  for proteomics and lipidomics alike; its replay-identity expectations pass on base. The other base failures are R1's.
- Head (fix): every test passes: 60 expectations in the new file, 312 in `test-prot-ruv-factors-in-model.R`, 1,153 in
  `test-omic-norm-ruv-design.R`; `compare.tsv` verdict `fixed` for all three files. A1 is met for this commit; R2 repeats
  it on the merged base.

## Decisions and deviations

- **Marker location (coordinator's ruling; ticket r1.4 7.1).** The DA is not a record step; the marker is the DA model
  record's `treat_moderation = "as_main_test"`, written while the current rule is in force by
  `acceptWorkflowOmicsDaManifest` and `omicsDaOperation` (lipidomics, metabolomics; legacy and declared lipid/metab
  designs) and by the proteomics legacy route's `da_model` (`differentialAbundanceAnalysisHelper`, `ruv_model` NULL).
  The DA replay (`replayWorkflowAnalysisRecord`) sets `MultiScholaR.treat_moderation_rule` from it: absent means the
  earlier rule (`"constant"`).
- **Declared proteomics designs and DIA with limpa** record no marker. Their computation does not read the option
  (`runProtRuvModelContrasts` passes the settings itself; limpa runs no TREAT), so their replay is unchanged
  whatever the option holds (T-A4 (c), (d) assert it under both rules).
- **Not rewritten:** `test-workflow-lipid-declarations.R:49` and `test-workflow-metab-declarations.R:49` (no
  declaration's `fixed` changes); `test-prot-ruv-factors-in-model.R:332-333` (both sides of that identity are
  `runTestsContrasts`, so it holds on head unchanged).
- **Rewritten** (6.2 (a)(ii)): `test-prot-ruv-factors-in-model.R` legacy digests (now: the legacy DA at 0 equals
  `eBayes(trend, robust)`, to 1e-10, for entry, entry with technical replicates and kernel) and its `da_model`
  expectation (gains the marker); `test-omic-norm-ruv-design.R:17` (0.385 to 0.389).
- **Continued session (T-A5 (c)).** The operation identity carries the build (`context$software`), not the rule; in the
  application a new build never serves an earlier build's product. Under `load_all` the build string is equal across
  commits, so T-A5 (c) models the new build by changing `context$software` before the re-run.
- **KIT-190 at merge.** Both branches append a row to `workflowReplayRuleMarkers`: a textual conflict at the list's
  closing line; keep both rows. KIT-190's `completeWorkflowAction` change filters its own marker from a declaration's
  `fixed`; this ticket's marker is in no declaration's `fixed`, so the two do not interact.
- **Old-record fixtures.** One committed record of a proteomics legacy DA made before the change exists:
  `tests/testdata/go-propagation/replay-base/record.json` (DIA, `DA-contrast-output`, cutoff 0; its model has no
  `treat_moderation`). Its project store is not committed; KIT-173's T7(b) replays it only inside that ticket's script
  `t7-replay.R` (skipped in the suite). Under this change it would replay with the earlier rule, so its DA, and the
  enrichment built on it, are expected unchanged (inferred; not run here). No committed lipid/metab record with a cutoff
  above 0 exists. T-A5 builds its records in process under each rule.

## Moved stored references

[R3's base-against-head list, against 6.2's prediction, to be filled after R3.]

## Figures measured on the earlier TREAT moderation (D9; labelled, none re-measured)

- KIT-124 S1 arms 1 and 3 and S1's stress cells (legacy proteomics test, TREAT at 0, one prior variance).
- KIT-177 R0 (run 45222712083; `kit-177/r0/effect.tsv`) and its Methods limitation.
- VAL-108's simulation (KIT-187 S2).
- KIT-185's legacy proteomics DA.
- KIT-197 (the pattern it explains was seen on the earlier test).
- The RUV panel's classic Legacy row before this ticket (KIT-124 arm 3: 0.385, 9.8 %); restated here from arm 28.

## For KIT-188: Skip against declared `none`

[Not yet determined: needs a run on the head after the merges (R2's commit). Recorded here when measured.]

## Not checked

- No R was run on this machine; every result above is a runner's.
- R0, M1, R2, R3: not run (wait for the KIT-196 and KIT-200 merges).
- The managed worker path's handling of the marker (R2's scope covers the worker tests).
- Windows; real data.
