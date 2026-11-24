# KIT-177 outcome: plain log2 at the label-free and TMT protein design step

Done and ready for its implementation audit (final section). The first part of this file is the record as it stood at
the F3 stop; the owner's decisions VAL-79 and VAL-80 then let the work continue.

Ticket revision 1.6 (step 0 applied 2026-11-24 in the main checkout). Branch `kit-177-prot-plain-log2` from
`v0.9-staging` fc48b6ec (KIT-180 merged). Runs on `run-tests` (R 4.6.1, reference BLAS, ubuntu-22.04); ids, commits and
job minutes in `runs.tsv`.

(Superseded; the record at the F3 stop.)

| Commit | Content |
|---|---|
| 025c9542 | R0 scripts only (`route.R`, `effect.R`, `replay-check.R`, `gui-entry-control.R`); `R/` and `tests/` equal to fc48b6ec |
| (records) | this file, `runs.tsv`, `feasibility.md`, `r0/` |
| (not run) | the tests commit and the fix commit, prepared and pushed but never sent to a runner (see "State of the branch") |

## Observed

(Superseded in part; the record at the F3 stop.) **R0 (run 45222712083, 122.9 job-minutes = 8.19 core-hours; expected 8.9, ceiling 18).** Three scripts, all `ok`.

- **F4** (S2 on `prot_lfq`, head 025c9542 against base fc48b6ec): the mechanism works from a script; every comparison
  empty, as expected with equal code (`feasibility.md`, `r0/f4-*.tsv`). R4f is not needed.
- **F3** (`r0/replay.tsv`): on the MaxQuant (LFQ) and Proteome Discoverer (TMT) e2e fixtures the records replay in
  order (history equal, the DA re-run), 5 of 6 later states `identical()`; `correlation_filtered` differs in
  `@args$daAnalysisUI` only; the DA long table is not `identical()` (its column names differ between the two reads);
  the kit fingerprint finds no mismatch. **Outcome (3), "replays in order, not identical"**: the ticket's stop rule
  (revision 1.6, from AUDIT-KIT-180-impl I06) stops here for the owner.
- **9.1** (`r0/effect.tsv`, per run `r0/effect-runs.tsv`, pairs `r0/effect-pairs.tsv`, fixtures
  `r0/effect-fixtures.tsv`): 400 of 400 route runs ok; first seed 4.92 minutes, projection 103.4 of 187; done in 103
  minutes. Means over 20 seeds (range in `effect.tsv`); bias = mean of sign(truth) x estimate over changed proteins
  tested (1 = unbiased); FDP and true-positive share at BH q < 0.05 on the reported p; "low" = lowest third of true
  abundance. `old` = log2(x + 1), `new` = plain log2.

| m | RUV | bias old / new | bias, low third old / new | FDP old / new | FDP, low third old / new | TP share old / new | TP share, low third old / new | calls differing |
|---|---|---|---|---|---|---|---|---|
| 10 | automatic | 0.815 / 0.941 | 0.551 / 0.745 | 0.190 / 0.174 | 0.646 / 0.711 | 0.528 / 0.456 | 0.145 / 0.084 | 42.5 |
| 10 | skip | 0.813 / 0.941 | 0.548 / 0.750 | 0.051 / 0.041 | 0.495 / 0.418 | 0.392 / 0.357 | 0.052 / 0.031 | 9.2 |
| 40 | automatic | 0.960 / 1.002 | 0.928 / 1.024 | 0.109 / 0.106 | 0.217 / 0.221 | 0.917 / 0.902 | 0.744 / 0.705 | 18.6 |
| 40 | skip | 0.956 / 1.002 | 0.919 / 1.023 | 0.043 / 0.039 | 0.115 / 0.107 | 0.895 / 0.881 | 0.684 / 0.647 | 3.9 |
| 150 | automatic | 0.984 / 0.993 | 0.980 / 1.000 | 0.113 / 0.112 | 0.148 / 0.151 | 0.969 / 0.968 | 0.901 / 0.898 | 4.3 |
| 150 | skip | 0.983 / 0.994 | 0.978 / 1.001 | 0.048 / 0.048 | 0.078 / 0.083 | 0.983 / 0.983 | 0.947 / 0.947 | 0.9 |
| 1,000 | automatic | 0.991 / 0.992 | 0.992 / 0.995 | 0.106 / 0.106 | 0.114 / 0.114 | 0.961 / 0.960 | 0.876 / 0.875 | 1.95 |
| 1,000 | skip | 0.991 / 0.992 | 0.992 / 0.995 | 0.049 / 0.050 | 0.050 / 0.051 | 0.979 / 0.979 | 0.934 / 0.934 | 0.05 |
| 10,000 | automatic | 0.992 / 0.992 | 0.994 / 0.994 | 0.110 / 0.111 | 0.117 / 0.117 | 0.961 / 0.961 | 0.876 / 0.876 | 0.75 |
| 10,000 | skip | 0.992 / 0.992 | 0.994 / 0.995 | 0.051 / 0.051 | 0.052 / 0.052 | 0.979 / 0.979 | 0.933 / 0.933 | 0 |

  The five e2e fixtures (both configurations): no call changes; largest change in a log fold change 1.5e-5 to 1.0e-4;
  on the two five-protein fixtures RUV automatic stops ("No valid percentage found") with either transform.

Inferred, not tested: at m = 10 and 40 the offset shrinks fold changes toward zero (bias 0.81 and 0.96) and plain log2
removes it (0.94 and 1.00); at low abundance plain log2 calls fewer true positives and, with RUV automatic, more false
ones in the lowest third, consistent with the unmodelled variance at low signal (ticket section 5, limitation). From
m = 150 on the two transforms agree to within the seed range. With RUV automatic the FDP is about 0.11 at every scale
for both transforms (the legacy RUV route's own behaviour, not this ticket's).

## Acceptance criteria at the F3 stop (superseded by the final table below)

| # | Status | Evidence |
|---|---|---|
| A1 | F1, F2 answered by reading; F3 answered (outcome 3, stop); F4 answered (yes) | `feasibility.md`, `r0/` |
| A2 | not reached (stopped before R1) | |
| A3 | not reached | |
| A5 | not reached (NEWS, Methods and ledger rows not written) | |
| A6 | not reached | |
| A7 | 9.1's table delivered (above, `r0/effect.tsv`); to be shown to the owner; T1b not run | `r0/effect.tsv` |
| A8 | not reached (R4 not run) | |
| A9 | not reached | |

(Superseded; the record at the F3 stop.)

## State of the branch

Prepared on the branch and not run (no runner has seen them): the tests commit (`tests/testthat/test-prot-design-log2.R`
with T1a, T2, T3, T3b, T4, T6, T7, T8, T9; the two contract rewrites in `test-workflow-native-session.R`; T5 not
written, because it depends on F3), then the fix commit (`R/`: the design step's rule, the counts and the refusal, the
notifications, the reopen notice, the `scale` warning, the texts) and the T1b script `t1b.R`. They wait for the owner's
answer on F3.

(Superseded; the record at the F3 stop.)

## Not checked

- Whether F3's non-identity is the moment of capture (instrument) or a replay difference (result): not tested.
- Why the DA long table's column names differ between the first run's read and the replay's read.
- The prepared tests and fix: never run; no claim about them.
- Proteome Discoverer's default "Reporter Abundance Based On" setting (F2: not found in Thermo's pages read).

## F3 diagnostic run (owner's decision VAL-79, 2026-11-24), declared before it is sent

- **Purpose.** Tell the moment of capture from a replay difference. The replayed `correlation_filtered` is read after
  the replayed DA (the replay store reopened at its end), as the first run's was read after its DA; the DA long table is
  compared column by column by name (rows aligned by protein), with both reads' column orders recorded.
- **Binding.** Run branch `run-k177-f3d` = `v0.9-staging` 55710a87 (R/ as shipped, KIT-172 and KIT-173 merged) plus
  the changed `replay-check.R`; nothing else differs. Script entry `timeout_minutes` 45; ceiling 1 core-hour (one
  scripts job: set-up 10 + about 5 minutes).
- **Decision (VAL-79).** If on both lanes the only differences are `@args$daAnalysisUI` (its timestamp) and the DA
  table's column order: T5 compares data slots and DA columns by name, the ticket says so (History, VAL-79), and the
  work continues. Any other difference: stop and report to the owner. One mechanics rerun at most (VAL-54); a
  comparison produced is a result and is not re-run.

**F3 diagnostic, run 45265502101** (commit 7194105d on `run-k177-f3d` = staging 55710a87 plus `replay-check.R` and
`route.R` of this branch; 11.7 job-minutes = 0.78 core-hours; `f3d/replay.tsv`). States: on both lanes 5 of 6 later
states identical; `correlation_filtered` differs only in `@args$daAnalysisUI`, and inside it only in `timestamp`, read
after the replayed DA as before it. The kit fingerprint: no mismatch. DA table: **not compared**: the replay's table
read back from the reopened replay store returned no rows and no columns (`da_rows 397 0`), so no column could be
compared by name. Diagnosis from the code: the script read the manifest with `hydrateSavedWorkflowEntity(...)$value`,
whereas the record reads it with `resolveWorkflowStateValue` (`R/utils_workflow_analysis_record.R:494-495`); on the
reopened owner the contrast entry it handed `readDIAContrastOutput` was not a contrast reference, so the long table came
back empty. This is the mechanics of the script's read, with no DA comparison produced.

**Mechanics rerun (the one VAL-79 and VAL-54 allow), declared before sending.** The only change: `readDa` reads the
manifest as the record does (`resolveWorkflowStateValue`), and the first run's DA table is read from its own store the
same way (the in-memory read recorded beside it), so both sides are read by one path. No expected value, comparison or
stop rule changes; the state comparison is not altered. Same binding form (run branch at staging 55710a87 plus the
script), same ceiling (1 core-hour), same decision as above.

**F3 diagnostic, mechanics rerun, run 45267978186** (commit 016b426a; 12.5 job-minutes = 0.83 core-hours;
`f3d2/replay.tsv`). On both lanes: replayed history equal to the record's; 5 of 6 later states identical, and
`correlation_filtered` differs only in `@args$daAnalysisUI$timestamp`; the DA long table read the same way on both
sides is `identical()` (397 rows, the same columns in the same order); the kit fingerprint finds no mismatch.
**Result: only the stored DA settings' timestamp differs.** Under VAL-79 the work continues: T5 compares the states'
data slots and `@args` without that timestamp, and the DA table by column name (rows aligned by protein); the ticket
says so (History, revision 1.6, VAL-79). This describes in-process records; nothing here is claimed for GUI-made
records (KIT-180 R0d).

## R1 declared (tests failing first)

Binding: `run-k177-r1` 8f7b4c6e = v0.9-staging 98df468b (differs from 55710a87 only in the methods ledger; R/ as shipped) plus the two tests-only commits (2cb0a7f1,
a9822043; nothing under R/). Files: `test-prot-design-log2.R`, `test-workflow-native-session.R`,
`test-workflow-analysis-record.R`; 3 shards, no base. Ceiling 4 core-hours (expected 1). Decides: T1a, T2, T3, T3b, T4,
T8, T9 and the two contract rewrites fail at the stated values; T5 and T7 pass (regression); T6 fails on its cell
rule (base: no change). A failure for another reason is recorded and, if it is the test's mechanics, corrected under
VAL-54 before the fix's run.

**R1, run 45270608493** (`r1/tests.tsv`, `r1/failures.txt`, `r1/kit177-lines.txt`; 3 jobs, 21.0 job-minutes longest,
48.1 job-minutes in all, 3.21 core-hours). As stated: T1a fails at (1), (2) and (3) for every shape, configuration and factor (e.g. TMT
`none`, s = 1000: design-step shift off by 1.73); T2 reports 2.88; T3 (0 not missing, NaN kept, log2(200 + 1)) and T4
(`log2(x + 1)`, no counts) fail; the contract rewrites fail (`log2(x + 1)`; 9.9955 against 9.9941); T6 fails its cell
rule (base: no change); T8 and T9 stop on the absent text function (no notification and no warning exist on base); T5
passes (regression). Two test-mechanics defects (VAL-54, instrument): T3b stopped because its fixture-session helper
was defined after the test; T7 failed only on matrix row names (the kept and design matrices are indexed differently).
Both corrected (tests only) and re-run as R1b on base.

## R2 declared (fix: targeted, dependents, R CMD check; T1b)

Binding: head c8fc140e (branch, fix and tests, staging 0.9.5 merged in) against base 55710a87 (staging), 4 shards per
side, `check` true; files: the three KIT-177 files and the ticket's regression scope with the dependents found by a
search for `buildProtDesignStateCheckpoint`, `protein_s4_initial`, `design_builder` and the changed functions
(26 files). Scripts job on the head: `t1b.R` (`timeout_minutes` 210; with `sva`). Ceilings: targeted 24 core-hours,
T1b 14. Decides: every new test passes on head; every other difference is a contract rewrite, a GUI-entry capture of
the three re-frozen lanes (until R4), or a stored reference listed for the re-record; check findings equal to base;
T1b's rule applied as fixed in the ticket.

**R1b, run 45274581502** (`r1b/`): T3b now fails at its stated value (the Spectronaut and stored-table zeros stay 0;
no count in the record); T7 passes on base, as a regression; every other result as in R1. Failing-first complete.

**R2, run 45274649782** (head c8fc140e against base 55710a87; `r2/compare.tsv`, `r2/check-compare.txt`, `t6.tsv`,
`r2/t1b.tsv`, `r2/t1b-seeds.tsv`). Tests: every new test passes on head (`test-prot-design-log2.R` 226 passed, 0
failed; T5 in `test-workflow-analysis-record.R` 351/0, count-changed by its 44 expectations); 23 files `same` (22 plus `test-workflow-native-session.R`) (base
failures in `test-prot-norm-core.R` 272 and `test-workflow-kit-harness.R` 1 identical on both); the contract rewrites
pass (`test-workflow-native-session.R` same, 734/0); `test-gui-entry-proteomics.R` 10 -> 64 failing: the 10 of base
(`prot_dia` 4, `prot_dia_limpa` 6) unchanged, and the three lanes `prot_tmt`, `prot_lfq`, `prot_lfq_fragpipe` differ
at every step from `design` on (18 steps each), none before it: the change the re-freeze (R4) absorbs. T1a: no
`KIT177-T1A-TP` line, so t and p held their tolerances in every configuration and factor. T6 (`t6.tsv`): on each
fixture the largest change equals its bound log2(1 + 1/min v); no call changes; largest fold-change change 8.5e-5 to
3.7e-4. R CMD check: 0 ERROR, 3 WARNING, 3 NOTE on both; none only on head.

**T1b** (scripts job, 96.5 minutes): `scale`: 14 of 20 seeds change a call (largest change in t 6.15) (in a fold change 1.27 log2 units) -> its own
ticket under the rule (D6); legacy RUV `automatic` 20 of 20 (t 10.5) (19 with a call change; fold changes moved by up to 971.6 log2 units, more than 5 in 46 of 80 seed and factor rows) and `manual` 20 of 20 (t 8.9) (fold changes up to 62.4 log2 units, more than 5 in 57 of 80) -> their own
ticket under KIT-172 D8's rule (compare VAL-81, KIT-188 widened to proteomics legacy RUV Automatic). TMT and LFQ shapes
give the same table. The new designs' rows did not run: every one stopped at s = 1 with "operation normalization-ruv is
missing parameter ruv_k_seed", a script mechanics defect (the tab sends the seed). Instrument run 1 of 3 (VAL-54): the
two new-design configurations only, with `ruv_k_seed` set to the package default (`RUV_FACTOR_COUNT_SEED`); no rule
or comparison changes; declared here before sending (binding the same head, ceiling the T1b row's remainder, 7.6
core-hours).

## R4 declared (control and freeze, S2)

Binding: S2 `gui-entry-control.R` on head e51eded6 (the fix; `R/` as c8fc140e) against base 55710a87 (staging,
containing KIT-180's re-frozen files), lanes `prot_tmt`, `prot_lfq`, `prot_lfq_fragpipe`, child processes under C
collation; `out` = `docs/qa/workflow-execution/kit/kit-177/control`, `out_pattern` `\.(json|tsv)$`, `timeout_minutes`
140. Ceiling 8.7 core-hours. Decides (section 6): step 1 frozen against base: no row on any lane; step 2: no row before
`design`; step 3: every row from `design` on classed `plain_log2 (KIT-177)` (values in states, artifacts and outputs on
keys both sides hold), none `unexplained`; state tables: at `protein_s4_initial` every cell base minus head =
log2(1 + 2^-head) to 1e-12, no other missing-value change; later tables identical in shape and identifiers. Only then
the three oracles are committed (scripts.tsv `ok`, no `withheld.txt`). Anything else stops the ticket.

**T1b, new designs, run 45290294513** (`t1b2/`; combined table `t1b.tsv`): `groups_all_features` and
`first_pass_nulls` (gate on): 0 of 20 seeds change a call, k or the control share at any factor; largest change in t
1.8e-11, in logFC 2.8e-14. With the legacy rows of R2, `t1b.tsv` is T1b's table: `scale` and both legacy RUV modes
reach "own ticket" by the fixed rules; the new designs are unit-invariant on this fixture.

**R4, run 45290288889** (head 32ead4fc, `R/` as c8fc140e, against base 55710a87; `gui-entry-control.tsv`,
`control/verdict.tsv`, `control/gui-entry-steps.tsv`, `control/gui-entry-tables.tsv`; 21.0 job-minutes). Both child
processes under C collation; `scripts.tsv` `ok`; no `withheld.txt`. Step 1 (frozen against base, KIT-180's re-frozen
files): 0 rows on each lane. Step 2: 0 rows before `design`. Step 3: from `design` on, 96 (`prot_tmt`), 102
(`prot_lfq`) and 102 (`prot_lfq_fragpipe`) rows, all `plain_log2 (KIT-177)` (288 `states` digests, 12 `artifacts`),
0 `unexplained`; no output, tab status, history or notification differs. State tables, cell by cell: at
`protein_s4_initial` every cell base minus head equals log2(1 + 1/v) to 1e-12 (0 violations, no zero turned missing,
no other missing-value change); later tables of the same shape and identifiers, largest absolute difference 2.4e-4 to
4.8e-4. The head record equals the file frozen from it. **Control passed**: the three oracles (source SHA 32ead4fc,
written under C collation by `guiEntryExpectLane`) are committed to `tests/testdata/gui-entry/` (SHA-256 prefixes
`prot_tmt` f35ae4853f91, `prot_lfq` 95553da3647a, `prot_lfq_fragpipe` dd1ed8e188d7, equal to the artifact's). Chain:
this record follows KIT-180's (`kit-180/gui-entry-control.tsv`).

## R3 declared (full suite, base against head, with check)

Binding: head c06568f3 (the fix and the re-frozen files) against base 55710a87, `tests` all, 20 shards per side,
`sva` on both, `check` true, `test-gui-entry-proteomics.R` 120 minutes. Expected 62.3 core-hours, ceiling 72. Decides
A3: every difference is a new test, a contract rewrite named in section 6, the re-frozen GUI-entry references (expected
`fixed` on the three lanes), or a TMT/LFQ stored reference listed for the re-record with its size; check findings equal.


**R3, run 45303214444** (head d7a04a08 against base 55710a87, 20 shards per side, `sva` on both; `r3/compare.tsv`,
`check-compare.txt`; 987.4 job-minutes = 65.83 core-hours, expected 62.3, ceiling 72). 454 files: 452 `same`, among
them `test-gui-entry-proteomics.R` 183 passed / 10 failed on both sides (the three re-frozen lanes pass on head; the
10 are `prot_dia` 4 and `prot_dia_limpa` 6, unchanged), the contract rewrites in `test-workflow-native-session.R`, and
KIT-172's four files; 1 `only-on-head` (`test-prot-design-log2.R`, 226/0); 1 `count-changed` (T5 in
`test-workflow-analysis-record.R`, 351/0); 0 `new-failure`, 0 `message-changed`. R CMD check: 0 ERROR, 3 WARNING, 3 NOTE
on both; none only on head. No stored reference outside the re-frozen files moved, so there is no re-record row for plan
section 11 beyond the ticket's own (the two contract rewrites; the three GUI-entry files re-frozen here, D2).

## Final: acceptance criteria

| # | Status | Evidence |
|---|---|---|
| A1 | met: F1 and F2 by reading, F3 outcome (3) then the VAL-79 diagnostic (only the DA settings' timestamp differs), F4 yes | `feasibility.md`, `r0/`, `f3d/`, `f3d2/` |
| A2 | met: T1a, T2, T3, T3b, T4, T8, T9 and the contract rewrites fail on base (R1, R1b) and pass on head (R2, R3); T5 and T7's steps pass on both; T6's bounds delivered; no T1a exceedance in t or p | `r1/`, `r1b/`, `r2/compare.tsv`, `r3/compare.tsv`, `t6.tsv` |
| A3 | met: targeted and dependents (R2) and full suite (R3) base against head: every difference is a new test, T5, or (R2 only, before the re-freeze) the three lanes from `design` on; check findings identical | `r2/`, `r3/`, `check-compare.txt` |
| A5 | met: NEWS under the unnumbered heading with 9.1's numbers (VAL-80) and F2's "some exports"; Methods sentences in `wording.md`; no version change | `NEWS.md`, `wording.md` |
| A6 | met: T8 (the design tab's start and the notifier the load-session handler calls; raised once for `log2(x + 1)` and for no transform; none for a head, DIA or other string); no output id added (R4 step 2 and 3: no `outputs` key on one side only) | R2, R3, `gui-entry-control.tsv` |
| A7 | met: 9.1's table shown and accepted (VAL-80); T1b's table with the rules applied | `r0/effect.tsv`, `t1b.tsv` |
| A8 | met: R4 passed before the three files were committed | `gui-entry-control.tsv`, `control/` |
| A9 | met: T9 (in process and through the finishing function the managed path uses): one warning and one log line with the count when `scale` runs on log2 values below 0; none for cyclic loess or for `prot_lfq_maxquant`; the result identical to a run without the warning code | R2, R3 |

**Core-hours** (`runs.tsv`): 109.87 used: R0 8.19, F3 diagnostic 0.78 and its mechanics rerun 0.83, R1 3.21, R1b 1.55,
R2 21.85, T1b instrument run 6.23, R4 1.40, R3 65.83. Expected 100.2 (revision 1.6), ceiling 245.4; every row inside its ceiling except R1, whose row (ceiling 4) holds R1 3.21 and its VAL-54 instrument run R1b 1.55,
4.76 in all, inside the ticket's total as VAL-54 (iii) requires (targeted 15.35 of 24; T1b 12.66 of 14; R4 1.40 of 8.7; R3 65.83 of 72) (main jobs only; with the plan and results jobs about 0.2 core-hours more).

**T1b's rules** (for the owner): `scale` changes a call in 14 of 20 seeds, so by D6's rule it is its own ticket
(removal or repair of the option); legacy RUV `automatic` (20 of 20) and `manual` (20 of 20) reach their own ticket by
KIT-172 D8's rule (VAL-81 widens KIT-188 to proteomics legacy RUV Automatic); the new designs are unit-invariant here.

**Gate-off and bundle lane.** The ticket defines no gate and no retained path (section 7: "No gate and no retained
path") and no bundle-lane run; neither was run.

## Not checked

- GUI-made label-free and TMT records' replay: not claimed (KIT-180 R0d); T5 and F3 are in-process records.
- The load-session site of outcome 4 is shown by T8 structurally (the handler calls the notifier on the restored
  owner); the handler itself was not driven with a saved session.
- T9's managed path is exercised through `finishProtNormNormalization`, the function the managed dispatch finishes
  through, not through a worker.
- Proteome Discoverer's default "Reporter Abundance Based On" setting (F2: not found in the pages read).
- Windows and the installed application: not run (no bundle lane in this ticket).

## Audit text repairs (AUDIT-KIT-177-impl-20261124, I01, I02, I04, I06, I07, I11), no code change

NEWS, the ledger row and `wording.md`: `scale`'s unit dependence on any data (I01; the in-app string is recorded in
`wording.md` as a known false sentence with the corrected proposal, for the owner's wording pass); T1b and T6 in the
ledger row with the fold-change sizes and the failed advance inference (I02); the low-abundance figures for both RUV
arms, the Benjamini-Hochberg sentence, "no call changed" (I04); records (I06); sessions adopted from earlier versions
not told (I07, NEWS); F1(a)'s scope completed (I11, `feasibility.md`). Commit messages: 9c8cb8ca added the ledger row,
a17c61d1 added NEWS (I12).

## R2b declared (on the merge with current staging; AUDIT-KIT-177-impl I05)

Binding: head = the merge 88892e21 (staging 17365cb0 merged into the branch; conflicts only in `NEWS.md` and the
ledger) against base 17365cb0; R2's 26 files plus KIT-187's `test-prot-missing-value-counts.R`; 4 shards per side;
`sva`; `check` true; `test-gui-entry-proteomics.R` 120 minutes. Ceiling: the ticket's R2 row, 24 core-hours (expected
16). Verdicts that pass: `test-prot-design-log2.R` `only-on-head`; `test-workflow-analysis-record.R` `count-changed` (T5);
`test-gui-entry-proteomics.R` `same` at 183/10 on both (the three re-frozen lanes pass on the merged tree with KIT-183,
KIT-186 and KIT-187 under them); `test-prot-missing-value-counts.R` `same`; every other file `same`; check findings
equal. Any other verdict is a result: stop and report.

**R2b, run 45313880809** (`r2b/`; 292.8 job-minutes = 19.52 core-hours, ceiling 24): every verdict as predicted. 25
`same`, among them `test-gui-entry-proteomics.R` 183/10 on both (the three re-frozen lanes pass on the merged tree with
KIT-183, KIT-186 and KIT-187 under them; the 10 are `prot_dia` and `prot_dia_limpa`) and `test-prot-missing-value-counts.R`
67/0 on both; `test-prot-design-log2.R` `only-on-head` 226/0; `test-workflow-analysis-record.R` `count-changed` 351/0.
R CMD check 0/3/3 on both, none only on head. The merged `R/`, `tests/`, `inst/`, `DESCRIPTION` and `NAMESPACE` are
staging's plus this branch's tested changes (the only hand-resolved conflicts were `NEWS.md` and the ledger). Core-hours
in all: 129.39 (109.87 + R2b 19.52) against the ceiling of 245.4.
