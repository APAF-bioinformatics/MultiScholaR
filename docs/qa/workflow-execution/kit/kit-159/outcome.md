# KIT-159 outcome (running record; build phase)

Ticket: `tickets/kit/KIT-159.md` revision 1.7 (1.6 plus the third re-check's one note, applied as step 0:
the label of 2.1's decision rule; no contract change). Branch `kit-159-validation-harness` from
`v0.9-staging` at a60c4e9a (version 0.9.1). Worktree `/home/doktersmol/Projects/APAF/MultiScholaR-kit-159`.

## Scope of the build phase

Steps that need only the code (section 11; the auditor's T0, T1, T2, R1): the harness `msrval` with its
unit and classification tests (tests first, shown failing on a runner: T0; then passing: T1, with
`covr`, `lintr` and the f_sim timing), the coverage regimes T-COV (T2, sized on the measured f_sim),
and the first record request (R1). Not in this phase: G1D (another agent, branch `run-k159-g1d`),
G1(a) at the release commit, G1(r), S0, L1a, L1b and later, P1 (they wait for KIT-172, KIT-176,
KIT-177, KIT-180 and KIT-182).

## Work items, in order (each: tests first, then implementation)

The harness is built in three module groups so that each group's tests are shown failing on a runner
before its implementation (no R runs on the laptop, so a single tests-only commit for every test of
section 10 would be one unreadable mass failure). T0 and T1 are each split across these groups; their
jobs are counted against T0's and T1's ceilings (6 core-hours each), not against new budget lines.

| Group | Files | Tests |
|---|---|---|
| a. scoring | `score.R`, `metrics.R`, `intervals.R`, `baseline.R` | 10.1 (all hand cases), T-VIEW, T-KEY, T-BLIND-1 |
| b. records and guard | `inputs.R`, `licence.R`, `declaration.R`, `provenance.R`, `determinism.R`, `evidence.R`, `units.R` | T-IN, T-LIC, T-LOG, T-DECL, T-PRE, T-ENV, T-SCORE, T-ANN, T-PROV, T-DET, T-CHUNK, T-RUN, T-RNG |
| c. driver | `driver.R`, `arms.R`, `parity.R`, `integration.R` (Part B later) | T-DRV, T-HELP, T-CHAIN, T-KER, T-TAB |

## Resume here (state 2026-11-24, builder A after the re-audit; builder B on `kit-159-replay`)

**Builder A (this branch; register `AUDIT-KIT-159-impl-recheck-20261124.md`, rulings VAL-137 and the
coordinator's two confirmations of 2026-11-24).** Builder B works on branch `kit-159-replay` (replay
of declared cells I07/E5, driving `val-run.R` E2, the G1(a) entry script and the D7 environment
recorder) and merges into this branch at its end; E4's `helpers_extracted` is B's (it is filled in
`val-run.R`). Shared REPAIR line: builder A 4 core-hours of it, plus DI10b's 1.5.

| Item | State |
|---|---|
| I10 record (VAL-137 (1)) | Corrected below and in `runs.tsv`: **by its declared rule the I10 diagnostic was "partial", not "confirmed"**. M1 and M2 held and every one of the 11 rows that failed in T2 rose (lowest 0.9605), but C4-145 `sens` (C4 g = 6) covers 0.7265 at every level: sensitivity is 1 in the full run and every deletion in 27.35 % of simulations (target 0.9992), so its interval has zero width there. T2's band (the declared band: fail below level - 0.02) reads that as a fail; the code labelled it `degenerate` by r1.13's rule, which the declaration did not name and T2's code did not have. Owner's ruling VAL-137 (1): the replacement of the 64 rows stands; C4 g = 6 `threshold_count` is recorded as not offered (no validation ticket uses it); the record states the partial outcome |
| E3 (`USAGE.md`) | done: consumer bullets for C4 and C4e to the re-pinned statuses; the "have not run end to end" line replaced (`valRunPhase`/`valScorePhase` ran end to end in R4d; `val-run.R`'s branches have not, E2); the duplicated loader sentence replaced by the built state; VAL-133 marked built; the stale "I10 not yet run" bullet removed |
| VAL-137 (2) family bound | code (`valWriteResults`), tests, `USAGE.md`, ticket 4.4: each component at its one-sided 95 % Clopper-Pearson upper bound on its miss rate, one convention |
| VAL-137 (3) no widening over independent units | code (`valLevelOne`, `valWidenFromRequests`), tests, `USAGE.md`, ticket 4.4: never widened at any count; each unit's within-run interval reported (also after a refusal fill); `results.md` states "over repeat runs, for the features scored"; widening only for jackknife forms (ii) and (iii) |
| E8 | `se_inflation` in the coverage rows, the summary and the pinned table (64 measured rows: half_width_ratio x z / t(g - 1), computed here by a standard-library t quantile checked against R's to 1e-12 and asserted against `qt` in R-E8); ticket's "1.00 for a correctly scaled interval" corrected |
| E7 | `generator` in the pinned table (`fixed_symmetric` for the 64 re-pinned rows, `redraw` elsewhere), in `valClassify`, the verdict table and `results.md` with its estimand |
| E4 | per-unit truths hashed into provenance (`truth.per_unit`, checked where recorded); `helpers_extracted`: left to builder B (`val-run.R`) |
| E6 | the licence guard's refusal tests assert the rule's message (R-I08, T-LIC-1 to 11, the hand-off and ticket-coverage cases, result tables, phase pre); no code change: the messages already name the rule; these tests pass on the old code by design (strengthening, not a fault demonstration) |
| DI10b | declared before the request in `t2/di10b-declaration.md` (what it measures, N1/N2 mechanism checks, which result means the symmetric split was needed); generator `fixed_onedir`; result below when in |
| Coordinator's confirmations | (1) under VAL-133 (2) the value path is taken only when the DA artifact is the only differing capture, any other difference `not_certified`; (2) a link's change classes default to the record's classes whose every part carries the ticket's tag: both in `USAGE.md` and the ticket |

**Failing-first at the audited fault, per finding (re-audit section 1; this replaces "every new test
fails" in the R0 and R11 lines below).** Shown failing AT the fault before the fix: I01 (R0 T-PH-5/7/8),
I07 (T-PROV-1 after its correction, R6), I15 (R0 R-I15; R5, R8, R9), I16 (R0), I18 (R0); as absence
of the mechanism: I03, I11. Shown only by the fix passing (the new test first failed at a
precondition or an API change, not at the fault): I02, I04, I05, I06, I08, I09, I13, I17 (and I20
partly: the anchor case at its fault). Of this builder's items: VAL-137 (2), (3), E4, E7, E8 and the
DI10b generator are shown in run A1 (below); E6 is a strengthening of passing tests.

**Runs of builder A** (REPAIR line unless stated):
- A1 (k159-rep-a1-d1c8d4e6, run 45334780977, tests only at d1c8d4e6): sent; purpose: show each new
  test failing at its fault; ceiling 25 min; decides nothing else.
- A2 (fix at 647706dd, same files with `--lint`): to send after A1.
- A2 (45336510695 at 647706dd): 49 of 50 pass, lint 0; R-V137-2 failed on a fixture defect (corrected in the test at 820f3207; re-run owed). 0.71 core-hours.
- DI10b (45336507597 at 647706dd, line DI10b; sent with A2): **valid; by its declared reading, conditioning suffices for C4e g = 5.** N1 (bias target -0.234) and N2 (rejection target 0.428) hold, so the one-directional shift is present; every statistic passes at every level (95 %: mae 0.9575, bias 0.960, rej 0.9485; no zero-width intervals). Compared with D (fixed_symmetric) coverage is lower but within the band (rej at 95 %: 0.9485 against 0.9835). Reported to the owner; it decides nothing by itself. 1.10 core-hours of 1.5. Tables `t2/di10b-*`.


Ticket r1.16 (main checkout; section 31). Rulings: VAL-109; KIT-199 plan audit ruling C1 (a); KIT-199
re-check Q01; KIT-192 implementation audit I04. Repair line: 15 core-hours (VAL-109).

**Observed (runs; tables in `rep/`, `runs.tsv`):**
- R0 (45310054283, 4b4e4355, tests only) and R1 (45310058037, 4a57cbdb, whole suite) were sent
  together, so the tests-first order rests on the commits (05d300b2/4b4e4355 add no `R/`), not on
  the send order. R0: every new test fails; most at the audit's fault (T-PH-5/7/8 I01; R-I03, I04,
  I06, I09/I18, I11, I15, I16, I20; T-TAB-OFF I09); some at a precondition (R-I01, I02, I05/07/08,
  I08: valPrepare; E2E: no da.tsv; T-PRE-RUV: `ruv` argument). Passing by construction: T-PROV-1's
  nested fields and R-I17's valRunCells backstop; both corrected (145621ac) and re-run on the
  tests-only code (R6, run branch `run-k159-r0fix`): T-PROV-1 fails at I07; R-I17's capture cannot
  be reached there (that valRunPhase cannot run a ticket of the repaired contract), so I17 is shown
  on the fix only.
- R1: 216 of 226 pass. All ten T-TAB routes, T-TAB-RUV, T-TAB-OFF on DIA_LIMPA and trx_screen (RUV
  off completes, before KIT-192) and T-PRE-RUV on DIA, LFQ, TMT, DIA_LIMPA pass on the merged
  product. Failures: E2E x3 and R-I01/I02/I05-08/I08/I06/I20 (fix defects, see below), T-SCORE-1
  (fixture mechanics), 27 lints in the repair code (a result; corrected cb3a57e9).
- R3 (4e78b246), R5 (96d87c12), R8 (20028e18), R9 (1d77f70a, pending): tests first for passes 2
  and 3 and the coordinator's items; each fails at its fault; R3's plan-generator test passed by
  construction, corrected (ee445cfc) and shown failing in R7 (run branch `run-k159-r3fix`). R3 also
  shows the I10 mechanism directly: under T2's generator the plain arm's null rejection share is
  0.37 (target 0.05).
- R4b (45315050609, 169aee9f): 50 of 55 pass (repair2, repair3, coverage, coverage-rank, chain,
  phases, most of repair); fails: E2E x3 (verdict row: row names from a named family vector),
  T-SCORE-1 (zero-row comparisons cbind), R-I08 (toy metric lacks `...`: mechanics), 2 lints.
  Corrected at 5901132a and 36297e25; R4c sent.

**Fix defects found by runs, each its item's one correction:** valPrepare `file.exists(NULL)` for
units without files (R5/R1, I02; 55e349a2); verdict row names (R4b, I05; 5901132a); phase score with
no comparison rows (R4b, I06; 5901132a); lints (R1, R4b).

**Then (all observed):** R4c (36297e25): repair and records pass, lint 0; E2E stopped because the
test passed no `pkg_root` (mechanics, corrected 89094e2a). R4d (89094e2a): the three end-to-end
cell-to-verdict tests pass on the merged product (prot_dia 56/56, metab_combined 59/59, trx_screen
56/56). I10 diagnostic (45316883316 at 169aee9f, T2's line): by its declared rule "partial", not
"confirmed" as first written here (corrected after the re-audit, S1: M1, M2 hold; all 11 rows
failing at 95 % in T2 rise, lowest now 0.9605; C4 g = 6 `sens` 0.7265 at every level, a fail under
the declared band, labelled `degenerate` by r1.13's rule); re-pin ratified by the owner (VAL-137
(1)); re-pinned
(d0492872): 64 C4/C4e accuracy and threshold rows replaced, C4 g = 6 threshold `degenerate`;
half-width ratios 1.4 to 5.3 (D) and 2.0 to 3.2 (consumer rank rows, H). R10 (d0492872): 48 tests
pass, lint 0.

**HANDED OVER AGAIN (second builder, context long) at the records commit after 80d17087.**

Done by the second builder (observed; tables `rep/r11-*`, `rep/r12-*`, `rep/r13-*`):
- Item 1 (I15 remaining G1(a)), `R/g1a.R` and `gates.R`: `valReadControlLink` reads a ticket's
  `gui-entry-control.tsv` in all four committed forms (a/b of KIT-177, 180, 199; runner/oracle and the
  comment line of KIT-172; `frozen_vs_base`/`frozen-base` -> base differences with their values;
  `base_vs_head`/`base-head0` -> before and after the changed step; KIT-172's `base-head1` and
  `head1-head0` not part of the link) and `declared-keys.tsv` (field prefix dropped; raw line breaks in
  cells read); link (i) checks each own declared row's stored and live values against `frozen` and
  `expected` (`valDeclaredValuesOk`); `valAdoptReferenceRows` compares (lane, step, class) and key
  sets in either column form, a truncated KIT-172 row as a subset; D7 / VAL-133 (2):
  `valCompareDaValues`, `valReadDaLink`, `valG1aDaMode` (exact only on the link's control CPU models
  when they agreed byte for byte, with the ten environment keys of KIT-199's corrected I04 text equal
  and the lockfile unchanged), `valG1aReleaseValues` (identical masked captures first; else values;
  certified with the difference printed, or not certified); runners `valG1aRunLane` (child process,
  inside `test_that()`, installed through the extracted helpers or under `load_all`, G1D's DA copy hook
  for `prot_dia_limpa`), `valG1aPart1`, `valG1aPart2` (chain from `git log --first-parent` of the oracle
  since 12d68ca7, each merge's `kit-<n>/gui-entry-control.tsv`; C51 label via `not_refrozen`).
  Fixtures: the real records in `harness/msrval/tests/testthat/g1a/` with `SOURCES.tsv` (source path,
  commit, SHA-256; KIT-199 at 6b3e615a, the others at staging 0b926295).
- Item 2 (VAL-133 (1)): `valCompareParityDa` compares values under SWD-083 after the digest
  (`identical`, `within_tolerance`, `differs`; no 10-digit decision; refuses without L1's table, keys,
  q and fold-change columns); the measured KIT-199 R4 rows (P00025 fdr and P00146 log2FC cross the 10th
  digit; P00189 log2FC 5.2e-11 relative; calls identical) fail the old rule and pass the new;
  `PARITY_INSTALLER_EXCLUDED` and `valInstallerAcceptanceRows` exclude the six `prot_dia_limpa` DA steps
  by name.
- Runs: R11 (tests only, 8b52be7d) every new test fails at the fault (absent functions; the parity
  test at its status); R12 (696a8e11) 25 of 28 with three fix defects and 4 lints, corrected once
  (80d17087); R13 (80d17087): all 28 tests pass (159 expectations), lintr 0.
- Interpretation to confirm (coordinator/owner): VAL-133 (2)'s value path is taken only when the DA
  artifact is the only differing capture; any other capture difference from the re-frozen file is
  `not_certified` (the G1(a) text's "expected identical"). The default change classes of a link are
  the record's classes whose every part carries the ticket's tag (a ticket does not list them
  separately in its record).
- Not yet run against a real lane: `valG1aRunLane`'s child (only a stubbed runner is tested); the
  first real use is G1(a)'s own job, or a one-lane smoke run (lipid_canonical, both modes, about 10
  job-minutes) before it. Not built: the job-environment recorder for D7 (KIT-199 `control.R`
  `writeEnvironment` is the template), the G1(a) entry script that chains part (1), part (2) and
  `gates.md`; L1's committed DA tables for G13 (the parity request needs `table`, `keys`, `q_col`,
  `fc_col` per route, and `valCheckParityL1` does not yet check a table's hash).

**Next, in this order (each tests first, shown failing at the fault, then the fix):**
1. Coordinator's re-audit list (register `AUDIT-KIT-159-impl-recheck-20261124.md`) and the owner's
   VAL-137 rulings: E3 (USAGE.md consumer bullets to the re-pinned statuses; the "have not run end to
   end" line; the duplicated loader sentence); E8 (`half_width_ratio` divides by z while the jackknife
   uses t(g - 1): add an SE-inflation column with the t factor removed; correct the ticket's "1.00 for
   a correctly scaled interval"); E7 (a `generator` column in the pinned table, the generator named in
   results.md: 500-feature rows conditional on one feature draw); VAL-137 (1): the re-pin STANDS, the
   C4 g = 6 sensitivity row recorded as not offered, outcome.md, runs.tsv and the ticket say the
   declared rule's outcome was "partial" and why; VAL-137 (2): the family-wise bound uses a one-sided
   95 % Clopper-Pearson upper bound on each component's miss rate (in place of the point estimate,
   and of 3/n at zero misses), results.md's family line and the ticket; VAL-137 (3): a run-level
   interval over independent units is never widened by a within-run interval at ANY count (each
   unit's within-run interval reported; results.md: "over repeat runs, for the features scored");
   widening only for jackknife forms (ii) and (iii): code, ticket text, USAGE.md; E4
   (`helpers_extracted` filled; per-unit truths hashed into provenance); E6 (refusal tests assert the
   message and which licence rule fired); E2 (drive `val-run.R` lifecycle, counted, score on the toy
   ticket in a temporary git repository); in outcome.md, per finding, whether a failing-first
   demonstration AT the audited fault exists (auditor: only I01, I07, I15, I16, I18, and I03, I11 as
   absence; I02, 04, 05, 06, 08, 09, 13, 17 shown only by the fix passing).
2. DI10b (owner and coordinator, budget line DI10b 1.5 core-hours from the unassigned total, 81.9 of
   276 used before it): one confirmatory cell with the fixed feature draw and ONE-directional changes;
   declare here BEFORE sending what it measures and what result would mean the symmetric split was
   needed; reported, decides nothing by itself.
3. I07 / E5 replay of declared cells: copy the cell's project store, prune it with
   `kitPrepareReplayProject` (`tools/profiling/workflow_kit.R:1936`, read and extracted from the
   package checkout, never edited) to the record's import and design states, then
   `replayWorkflowAnalysisRecord(record, prepared, inputs, paths, scratch)` and `valReplay`'s
   comparison; replace `replay = list(status = "not_built")` (`runfolder.R:304`). The package's own use:
   `tests/testthat/test-workflow-lipid-declarations.R:198-216`.
4. E1 LAST, with covr: the whole suite with `--lint` and T-TAB at ONE commit merged with current
   staging (KIT-175, 177, 187, 189, 191, 192, 199 in it); repeated once after KIT-182, 190, 194 merge.
   covr estimate 3 to 5 core-hours (the plain suite 25 min on 4 cores; instrumented 2 to 3 times).
Owed after KIT-202 merges (not now): S0's three metabolomics routes and `metab_lc` export-then-load
re-run (about 1 core-hour). Do not repeat R1 for later product merges until the coordinator says.

**Lines:** REPAIR 14.81 of 25 (R0 0.86, R1 2.21, R3 0.70, R4 0.20 cancelled, R5 0.90, R6
0.68, R7 0.57, R8 0.78, R4b 0.87, R9 0.81, R4c 0.88, R4d 2.27, R10 0.77, R11 0.76, R12 0.77, R13
0.78); T2 29.00 of 39 (DI10 2.07); DI10b 0 of 1.5 (new line, coordinator 2026-11-24).
VAL-54 instrument runs per item: R0 item 1 (R6), R3 item 1 (R7), E2E item 1 (R4d). Item 1 (G1(a)
remainder): one correction used (R12 -> R13).
Inferred, not run: the whole suite at the final head (R1 ran at 4a57cbdb).

## Builder B (branch `kit-159-replay`): replay, `val-run.R` driven, G1(a)'s job pieces, E4

Branch `kit-159-replay` from `kit-159-validation-harness` at b8683b3d, worktree
`MultiScholaR-kit-159-b`; merged with builder A's heads 67594694 (f4a3152d) and 69bb2e4f (this
merge). This section replaces builder B's `outcome-b.md`, and B's runs are rows RB0, RB1, SM0, SM1
of `runs.tsv` (formerly `runs-b.tsv`).

| Item | Register | What was built | Tests | State (observed) |
|---|---|---|---|---|
| 1. Replay of declared cells (G5) | I07, E5 | `R/replay.R`: the cell writes its analysis record inside its session (`writeWorkflowAnalysisRecord`, inside the owner's state scope, 1080f82a); after the session a copy of its store is pruned to the record's import and design states by the package's `kitPrepareReplayProject` (extracted from `tools/profiling/workflow_kit.R` with its closure, never edited), the record replayed there with `replayWorkflowAnalysisRecord` on fresh experiment paths, under the cell's thread variables and seed; the two stores compared by the package's own fingerprints (every state, as KIT-028's `kitReplayJourney`, and the head's DA result). Statuses `identical`, `differs` (differing states and `da_results` named), `unsupported` (class); `not_declared`, `not_run`. `replay = not_built` removed; `valCheckProvenance` refuses a cell without such a status and a `differs` without its differences | T-REP-1 to T-REP-4, T-PROV-2 | RB0: each failed at its fault. RB1: all pass; **T-REP-4, a real `metab_combined` cell replayed: `identical`** |
| 2. `val-run.R` driven itself | E2 | The cell and score phases read the committed declaration through `valCheckDeclaration` first (not the working tree); any `val_error` ends the job with one template line (class) and exit 3 | T-RUN-E2 | RB0: lifecycle, counted, score ran to completion on the toy ticket in a temporary repository; the refusals failed (exit 1, no line). RB1: 21 of 21 |
| 3. G1(a)'s job pieces | I15 remainder, VAL-120 (2), VAL-133 (1) | `R/g1a-entry.R`: `valJobEnvironment` (CPU model, image, the ten D7 keys, lockfile SHA-256; derived as KIT-199's `writeEnvironment`), `valLockfileUnchanged`, `valG1aDaInputs`, `valG1aGatesMd`, `valParityL1Tables` (L1's DA tables for G13 read by path and checked by SHA-256); the entry script `val-g1a.R` (environment, part (1), part (2), `gates.md`) | T-G1A-ENV-1, -ENV-2, -GATES, -ENTRY, T-PAR-5; the smoke run | RB1: all pass. SM1 (real lane): `lipid_canonical` part (1) `identical`, installed against `load_all`, each in its own R process; environment and `gates.md` written (`g1a-smoke/`) |
| 4. E4 | E4 | `helpers_extracted` filled by `valRunPhase` when `val-run.R` passes none (GUI-entry helpers for workflow cells, the replay kit where it replays). Per-unit truths: builder A's `truth.per_unit` kept at the merge, B's duplicate removed | in T-REP-4, T-PROV-2 | RB1: pass |

**Lint (RB1).** `replay.R` 186 and 200, indentation only; corrected after this merge as item 1's one
correction (layout, no behaviour change); shown by the merge run below.

**Coordinator's rulings on builder B's two points (ticket r1.18, section 33; `USAGE.md`).**
- (a) ACCEPTED: G5 compares the package's content fingerprints of every state and of the head's DA
  result between the cell's store and the replayed store, because the replay re-runs the DA without
  publishing tables (proteomics `publish = FALSE`; the omics DA writes no files), so no replayed table
  text exists. At least as strict as 7.3's masked digest then value. `identical` = every state
  fingerprint and the DA result's fingerprint agree; `differs` names what differs and carries no
  largest value difference.
- (b) NOT decided, with the owner: whether L1 may write its driver DA tables (keys, q and fold-change
  per route on generated inputs) for G13, given 5.2's rule that no feature-level table leaves a
  lifecycle job. Nothing more is built for it; the reader and its SHA-256 check (`valParityL1Tables`)
  stay. **Until it is ruled, G13 cannot take its last link** (the bundle driver's DA table against the
  L1 driver's on `run-tests`, digest then 10-digit fingerprints): there is no committed L1 table to
  read, so P1 can establish browser-against-driver parity on the bundle, but not bundle-against-L1.

**G1(r)'s `g1a-reference.tsv`: not built (its content is not fixed unambiguously).** The row form
(lane, step, class, fields, keys; reader `valReadReferenceRows`) and the source (run 44953187101's
`failures.txt`, SHA-256 1b918062…) are fixed, and re-freezing tickets have committed rows for eight
lanes (KIT-172: `lipid_canonical`, `metab_*`; KIT-180: `prot_tmt`, `prot_lfq`,
`prot_lfq_fragpipe`; KIT-199: `prot_dia_limpa`). To be decided before it is built:
1. the class of `prot_dia_limpa`'s reference steps from `da` on: r1.16 withdrew
   `limpa_da_last_digits`, so by G1(r)'s text the step "fits no class" and G1(r) stops, while
   KIT-199's committed rows class the same reference-run steps `declared_change (KIT-199)` (a
   ticket's change class for a difference measured at 12d68ca7, before that ticket) and give two
   rows per summary step (`declared_change (KIT-199)` and `integration_artifacts`). G1(r) adopts a
   committed file only if it equals G1(r)'s own rows, so which class G1(r)'s own row carries decides
   whether it stops;
2. who builds G1(r)'s own rows for all ten lanes from the 69 blocks, and how: the harness has no
   `failures.txt` reader, so either a reader is written and run on a runner, or the rows are read by
   hand and recorded as such;
3. `g1a-reference-source.txt` (the 69 GUI-entry blocks with the source's hash) is still to be
   committed before the artifact expires (about 2026-12-21; the file is also kept in
   `tickets/kit/evidence-run-44953187101/`).
`val-g1a.R` refuses part (2) without the table (T-G1A-ENTRY), so part (2) cannot run until it exists.

**Builder B's runs** (REPAIR line; rows in `runs.tsv`): RB0 0.91, RB1 1.13, SM0 0.01 (refused at
plan: `requested_by` over 300 characters; VAL-54 instrument defect, 1 of 3), SM1 0.86 (SM0 re-sent
unchanged but for that line). The merge run: below.

**Not checked (builder B).** Part (2) of G1(a) on a real lane (needs `g1a-reference.tsv`); part (1)
on the nine other lanes; replay on any route other than `metab_combined`; G13's last link (ruling
(b)); the whole suite and covr (last, on the coordinator's word).

## Earlier resume record (state 2026-11-24, pre-release repair written, UNRUN)

The pre-release implementation audit (`tickets/kit/audits/AUDIT-KIT-159-impl-prerelease-20261124.md`)
found the cell-to-verdict path unfit (I01 to I20). The owner answered every waiting point in methods
ledger VAL-109 (staging 701458be). This section is the state of the repair; the next section is the
state before it.

**Commits (branch `kit-159-validation-harness`), every one UNRUN:**

| Commit | What | Runs it needs |
|---|---|---|
| e1a607a1 | merge of `v0.9-staging` e16ddbae (I12; KIT-180, 172, 173, 186, 183, VAL-90); one conflict, `docs/methods-decision-ledger.md` (staging's VAL-63 to VAL-107 kept in order, the branch's `VAL-<merge>` row after the last VAL row) | in R1 |
| 05d300b2 | tests first: `test-repair.R`, `test-e2e.R`, `test-tab-repair.R`, `helper-repair.R`; `test-phases.R` and `test-records.R` to the new contracts | R0 (shown failing) |
| 4b4e4355 | tests first: T-PH-8 to the O8 ruling | R0 |
| 8f37e3e1 | the fix (I01 to I09, I11, I13 to I18, I20; O1 flag, O7 (b), O8, D19) | R1 (passing), R2 (covr) |
| (records) | S0 ledgers committed (I19 (1)); `outcome.md` and `USAGE.md` corrected (I19 (2), (3)); `runs.tsv`: T1cov and T1lint relabelled as results (O6 (b)), the REPAIR line (VAL-109) | none |

**Runs needed, on the REPAIR line (15 core-hours, VAL-109), estimates from S0 and the suite run
(a classic proteomics corpus route run about 230 s; the suite before the repair 1,126 s on 4 cores):**

| Run | Binds | Purpose and what it decides | Estimate |
|---|---|---|---|
| R0 | 4b4e4355 | the new and changed test files only (`--files repair,e2e,tab-repair,phases,records`), shown failing where the audit's faults are; decides nothing else | 15 to 25 min, about 1.5 core-hours |
| R1 | the records head (harness = 8f37e3e1) | the whole msrval suite with `--lint` on the merged product: the fixes' tests and the I12 re-tests (T-TAB on LFQ and TMT after KIT-180, lipid and metab after KIT-172, DIA_LIMPA and trx RUV off, the end-to-end tests of three omic families); a failure in the mechanics of a test is VAL-54, a failing expectation is a result | 80 to 100 min, about 6 core-hours (`timeout_minutes` 150) |
| R2 | as R1 | covr per file at the merged head (I19 (4)) | about 2 core-hours |
| D-I10 | a commit with the generator change (not written) | the I10 diagnostic, four 500-feature cells at 2,000 simulations, on T2's line | about 2 core-hours (T2) |

Repair line: about 9.5 of 15 core-hours if nothing is resent. A re-run after the last pre-release
product merge (KIT-177, 182, 187, 190, 192, 194; I12) is a further R1, about 6.

**Waiting / not built (with the register's reference):**
- I10: the diagnostic generator (fixed per-cell feature draw, changes split +1/-1) and its tests are
  NOT written; the ruling (VAL-109 (2)) asks for its confirm/refute rule to be declared before the
  send. Nothing about C4/C4e changes until it runs.
- O4: the half-width ratio for consumer-used `conservative` rows: not built. The zero-miss rule is
  written in `USAGE.md`.
- O5: printing both labels on the 117 T2 rows: not built (`USAGE.md` states the mix); needs r1.13's
  rule applied to T2's detail rows to produce the second label.
- O9: ticket text of 1.2, 7.2 and G1(b) to amend (ticket, not code).
- I07 replay of declared cells: not built (provenance says `not_built`); needs G5's knowledge of
  `replayWorkflowAnalysisRecord`'s inputs.
- I15 G1(a) part (1)/(2) runners and the live-difference classifier: still to be written, not only
  run; `limpa_da_last_digits` is out of the default classes until KIT-199 merges (D19, VAL-109).
- Ticket 5.7 needs `ledger-<route>-<arm>-<unit>.tsv` as the ledger name (I07).

**Rulings applied in code (VAL-109):** O7 (b): after a refused unit a verdict only where a regime is
registered for n - k units, else `not_evaluable` with the reason in `regime`. O8: jackknife widened
by the `full` unit's within-run interval (same view); independent units report each unit's
within-run interval in `flags`, no widening. O1: bands as run; results.md prints measured coverage,
the family bound (3/simulations for zero misses; replaced by VAL-137 (2)) and a FLAG where the worst-condition miss rate is at
least twice nominal; `fail` is not read on the miss rate. I16 as the register recommends.

**Divergences to record:** per-unit truths are built in the scoring process, not as the unit's
first cell inside the chunk (I20); `results` off proteomics stays NULL (documented, I14); within-run
component outcome across units: the units' common outcome, else none (I16; O8 ruled only the
widening).

## Earlier resume record (state 2026-11-24, S0 done, whole suite passes; before the audit)

**S0 done: all 15 route runs valid, every step `ok`.** Evidence `s0/s0r1/` (run 45277101451 at
a29a43c0: the 9 non-classic-proteomics cells) and `s0/s0i2/` (run 45288286009 at df71922f: the 6
classic proteomics cells), each with `routes.tsv`, `steps.tsv` (head after every step: G1(b)'s
expected states), `alternatives.tsv`, `modal.tsv`, `corpus.tsv` (generator, arguments, bytes,
SHA-256 of every generated file), `run.tsv`, `outputs.sha256`; the ledgers stay in the runs'
artifacts (1.7 MB). Shown to the owner; gates nothing (VAL-41).

Observed (synthetic inputs at default QC, shipped arm):
- Every route reaches a DA table: prot_dia 659 rows / 44 called; prot_dia_limpa 794 / 36; prot_lfq
  787 / 61; prot_lfq_fragpipe 787 / 34; prot_tmt 789 / 49; trx_screen 1,668 / 129; lipid_canonical
  366 / 29; metab_lc 293 / 38; metab_gc 292 / 56; metab_combined 585 / 94. Classic proteomics
  normalisation reaches `ruv_corrected`; no refusal anywhere (the ITSD diagnostic step is `ok`).
- Export-then-load: identical DA table (masked digest) to the accepted script on `prot_dia`,
  `lipid_canonical` and `metab_lc`.
- DIA_LIMPA without "Compute RUV-III Covariates": DA runs from `correlation_filtered`; the same
  794 keys and 36 calls, 0 calls changed, largest difference of log2fc, p or q 0.153 (digest differs).
- Two-species `prot_dia`: the organism modal is raised and confirmed at its rendered values
  (`selected_organism` 9606, `filter_to_organism` FALSE, `mixed_species_fasta` TRUE forced); DA
  table identical to `prot_dia` on the same report.

Inferred (not tested): with the filter off, the mixed-species box changes nothing downstream on
this input; the covariates button moves values but not calls here, on one synthetic input only.

**Instrument runs (VAL-54), S0 item, 2 of 3 used:** the first S0 send stopped 6 cells after DA in
the harness's table reader (`valPreRuvState` reads `normalized`, which a run with RUV applied never
commits: T-TAB had only run those routes with RUV off). i1 (8fea7419, 45285928330): lineage-walk fix
did not hold; the gate stopped S0 before any route, as designed. i2 (df71922f): the
`ruv_corrected` state's record read directly; new test T-TAB-RUV passes; the 6 cells rerun
(rerunning only invalid parts). The 9 cells of the first send were not rerun: the fix changes only
the pre-RUV read of classic proteomics routes, a path none of them takes.

**Coordinator's rulings of 2026-11-24 (engineering, recorded):**
1. The whole msrval suite with `--lint` once at the head, charged to T1. **Sent as run 45294023250
   at 788903bb: no result.** The script reached its 17-minute limit (exit 124) inside T-TAB's
   metabolomics routes; no test table, no lint. The limit was set so the job stayed inside T1's
   1.91, on an estimate of 1.55 taken from the fastest earlier whole-suite job (19.6 min); earlier
   whole-suite jobs ran 22 to 31.5 min, so the estimate was wrong. Cost 1.84 (T1). The coordinator
   then ruled it an instrument defect (VAL-54) and had it resent with a 35-minute limit, charged
   to S0's line (which carries the tests of its instrument fix df71922f). **Resend: run
   45298108528 at 4e76758d (harness code as 788903bb): 204 tests, 1,617 expectations, all pass (0
   failed, 0 errors, 0 skipped); lintr no lints; script 1,126 s; 1.93 core-hours.** Evidence
   `suite/suite-*`.
2. The coverage table's half-width ratio is NOT built now: the register's premise was wrong (the
   detail tables do not hold the simulated statistics), and it would need coverage cells re-run,
   which is evidence. Ticket 4.4 keeps the text, T2's rows `not_recorded`. **Candidate for the
   harness auditor and the owner:** keep each simulation's half-width and statistic in
   `valCoverageRegime`, print the column, and decide whether any T-COV cell (first C4e g = 14,
   offered `conservative` at coverage 1.000) is re-run for it.

**Core-hours by line:** T0 5.95 of 6; T1 20.93 of 21 (19.09 + 1.84, the suite run); T2 26.93 of 54;
G1D 1.1 of 12; R1 2.78 of 4; **S0 9.65 of 16** (expected 6) = S0t0 0.70 + S0 3.55 + i1 0.92 + i2
2.55 + whole suite 1.93. G1R, L1a, L1b, L2, L3, P1, M1: nothing used.

**Bound to the release commit (not started; stop here):** G1(a) and its reference table G1(r)
(with G1R only if needed); `prot_dia_limpa` still awaits D19; L1a (proteomics lifecycle with the
two-species route, the G1(b) alternatives at the release commit, G6 third run); L1b (other omics,
guard, failures, chunks, G14 (b) and (c), phase `score`); L2 (determinism twin, G8); L3 (check,
base against head); the parity stage's first real run (P1, `distributable.yml` stage `parity`,
Linux and Windows, 7.9); the MOFA lane (Part B, M1) with KIT-170; R1's live incomplete-answer
refusal with KIT-175's merge.

**Open points for the harness's pre-release implementation audit:**
1. (Closed) The whole-suite regression check: run 45298108528, all 204 tests pass, no lints.
2. The half-width ratio column (candidate, above).
3. The fix of df71922f (`valProtRuvInput`: the classic proteomics pre-RUV matrix read from the
   `ruv_corrected` state's normalisation product) changes what `matrix_pre_ruv.tsv` holds on DIA,
   LFQ and TMT with RUV applied; it follows the product's own integration-input reader. T-TAB-RUV
   checks columns and that it differs from the corrected matrix, not its values.
4. Divergences from the ticket text: corpus writer in msrval, not `lifecycle/`; FragPipe from the
   e2e generator's writer, not `lifecycle/corpus-fragpipe.R`; DIA from the e2e DIA-NN writer, not
   `workflow_perf107_synthetic.R`; metabolomics generated inputs are MS-DIAL (`msdial`), the e2e
   fixtures `custom`; ticket 1.2's "no FragPipe generator exists" is out of date.
5. Generated input sizes are declared, not tuned; the annotation rule on generated inputs uses the
   corpus's own synthetic UniProt cache (min share 1), never the network (ruling 4).
6. (Closed 2026-11-24, audit I19 (1)) S0's 15 ledgers are committed under `s0/s0r1/ledgers/` and
   `s0/s0i2/ledgers/` from the artifacts of runs 45277101451 and 45288286009; both folders now
   pass `sha256sum -c outputs.sha256` (the script-1 files `tests-harness-tests.tsv` and
   `tests-lint.txt` sit beside them, outside that list, as before).
7. The 9 cells of the first S0 send ran at a29a43c0, the 6 rerun cells at df71922f (the fix does
   not reach the 9 cells' paths: inferred from the code, not re-run).
8. The register's budget misreading (T0 "5.95 left" for 5.95 used) and the 4.4/10.2 placement are
   recorded in ticket section 29.

### S0 build record (2026-11-24)

**Ticket r1.14** (main checkout; register `AUDIT-KIT-159-build-points3-20261124.md` applied as
written, with two corrections recorded in its section 29: the failing-first run is charged to S0's
line, not T0's, by the coordinator's ruling, since T0 had 0.05 core-hours left; the half-width ratio
of point 4 sits in 4.4, where the "share of zero-width intervals" sentence is, and cannot be computed
from T2's detail tables, which hold per-row coverage, Monte Carlo SE and target (and, for T2R only,
the zero-width share; corrected 2026-11-24, audit I19 (2)) but no per-simulation half-widths or statistics: those rows print `not_recorded`; its code and any re-run
are the coordinator's decision, not built).

**S0's three pieces (r1.14), built tests first:**
- tests only 3e4b5e7a: `test-corpus.R` (T-CORP-1 to 6), `test-alternatives.R` (T-ALT-1 to 4,
  T-S0-1, 2); run 45275103241 (S0t0): all 12 fail as expected (missing functions; `alternative` an
  unused argument); 0.70 core-hours, charged to S0's line (coordinator's ruling). `s0/s0t0-*`.
- code a29a43c0: `R/corpus.R` (`valCorpusSpec`, `valCorpusFor`, `valWriteCorpus`, `valCorpusCell`:
  the five proteomics formats by the e2e fixture writers of `tests/testdata/e2e/generate_fixtures.R`,
  whose abundances are `tools/profiling/kit_protein_corpus.R`'s recipe, taken with
  `valExtractHelpers` so the file's top-level generation never runs; transcriptomics, lipidomics,
  metabolomics by `tools/profiling/kit_omics_corpus.R` through `sys.source` behind its own
  `sys.nframe()` guard; both read only, a test hashes them before and after; a synthetic UniProt
  cache from each FASTA as `tests/testdata/gui-entry/uniprot/build.R` builds it; the two-species input
  = the `prot_dia` corpus with the last 30 % of FASTA entries moved to E. coli K12, taxon 83333, and
  `mixed_species_fasta = TRUE` as a declared forced input); `scripts.R` (`valActionScript(...,
  alternative = "export_load" | "no_covariates")`); `R/s0.R` (`valCompareAlternative`, `valS0Plan`,
  `valS0TestsPassed`, `valRunS0`); `val-run.R s0`.
- **S0 sent at a29a43c0: run 45277101451** (request `k159-s0-a29a43c0`): script 1 the 12 tests with
  `--lint`; script 2 `val-run.R s0`, which refuses before any corpus or route unless every test of
  script 1 passed (VAL-54 instrument fault: fix and resend, at most three instrument runs). Script
  limits 20 + 190 min; with S0t0 the worst case is 15.4 of S0's 16 core-hours.

**Divergences recorded (for the auditor):** (1) the corpus writer lives in msrval (`R/corpus.R`),
not `lifecycle/` (7.2), so that it is tested, linted and pinned (conventions/R.md: a package, not
sourced scripts); (2) FragPipe is written by the e2e generator's existing FragPipe writer (KIT-101
r2.2, from the LFQ recipe), not a new `lifecycle/corpus-fragpipe.R`: 1.2's "no FragPipe generator
exists" was out of date; (3) DIA uses the e2e DIA-NN writer (KO/WT effect from the KIT-021 recipe,
UniProt-style FASTA that the organism modal reads), not `workflow_perf107_synthetic.R` (no group
effect; FASTA headers without OS=/OX=); (4) metabolomics generated inputs are MS-DIAL exports
(`import_tool = "msdial"`), the routes' e2e fixtures use `custom`; `metab_lc` runs LCMS_Pos only,
`metab_gc` GCMS only, `metab_combined` both.

**Budget:** T0 5.95 of 6; T1 19.09 of 21; T2 26.93 of 54; R1 2.78 of 4; S0 0.70 of 16 (S0t0) plus
run 45277101451 in flight.

---

Earlier state (R1 done; kept as written):

Ticket r1.13 (VAL-69). Every msrval test has passed (whole suite 45219635014, corrections
45225211841 and 45265533254; lintr 0). Rank re-run R (45227824000) done; USAGE and H2 updated
from it.

**R1 wiring (coordinator's word, 2026-11-24: charged to T1 as unit tests of harness code, VAL-69's
purpose; T0 was exhausted):** tests only 1992d121, run 45264824991 fails first (0.77); code
d6d3c1e1 (`valRecordUniprot`, `val-run.R record --uniprot`, `uniprot.tsv` on both allow-lists), run
45265533254 passes with lintr 0 (0.90).

**R1 done.** First send (45268141191, 45268145400 at 69b693af) stopped on a harness defect
(taxon passed as text; the product formats it with `%d`); by the coordinator's ruling an
instrument defect under VAL-54. Instrument run 1 of 3 at c08787c1 (fix `as.integer(taxon)`; no
expected value or rule changed): 45270635025 (`r1a`) and 45270640790 (`r1b`) completed. In both:
the open file at 1,086 bytes and SHA-256 566d7b3d... as r1.12; the UniProt record of the 25
accessions (taxon 9606) 25 rows, release 2026_04 read from `x-uniprot-release`, table SHA-256
d099e6e7...; library table byte-identical between the twins; phase pre `ref-qc-shape.tsv` (21 rows,
6 columns); provenance and `outputs.sha256` in every folder. Evidence `r1/r1a/`, `r1/r1b/`. The
twin comparison by `valCompareRecords` is L1b's step. The live incomplete-answer refusal waits for
KIT-175 (T-ANN-1 mocked stands).

**S0 not started (superseded by r1.14 above): preconditions do not hold.** S0 (10.4) needs the G1(b) alternatives: both QC
paths on label-free and TMT "the bypass still being present there", export-then-load, DIA_LIMPA
with and without the covariates button, and `prot_dia` on a two-species input through the organism
modal, all on generated inputs. On `v0.9-staging` the bypass is gone (KIT-180: no shell sets
`chain_rules$qc_bypass`), and the harness has no generated-input corpus writer, no alternative
scripts and no two-species input yet (`action-scripts.md` open points). Building them is new
harness code with tests first: T0 has 0.05 and T1 1.91 core-hours left. The coordinator decides
S0's form now that the bypass is removed, and its budget.

**Budget:** T0 5.95 of 6; T1 19.09 of 21; T2 26.93 of 54; R1 2.78 of 4.

**For the harness auditor (noted, nothing changed):** in the rank re-run, C4e g = 14 (rank,
delete-one-block) has simulated coverage 1.000 at every level (no miss in 2,000 simulations) with a
target of 0.049 (not at a bound, zero-width share 0), so it is `conservative` everywhere under
r1.13 and offered.

## Earlier resume record (state at dcc27fb5)

**Read first:** the brief `tickets/kit/audits/IMPLEMENTER-BRIEF.md`; the ticket `tickets/kit/KIT-159.md`
(revision 1.10, main checkout); "Rulings of 2026-11-23" below. Nothing runs on T1 and T2 stays held until
the owner answers the development-run rule; tests-only runs go under T0 (**3.60 of 6 core-hours used**,
after T0d).

**Built and verified on a runner:** group (a) scoring (T1a correction, run 45125695623).
**Built, verified except one code defect, fix committed unrun:** group (b) (T1b correction, run
45129042733; T-LOG: append mode, fixed in `valRunAsChild`).
**Built, fix committed unrun:** the coverage range split (T1s, run 45129513148: `spec[["target"]]`).
**Committed unrun:** group (c) code and the defect sweep (T0c, run 45131922735, showed its 25 tests
failing first); **group (d), the driver core** (T0d, run 45134811170, showed its 48 tests failing first;
code at dcc27fb5, below).

**Group (d) as committed (unrun), by file under `harness/msrval/R/`:**
- `scripts.R` (a file 7.2 does not list; 7.2 puts `valActionScript` in `driver.R`, split for
  readability): `valStep`, `valActionScript` (the accepted scripts of 7.4 for the ten routes; QC
  sub-tabs from the product's spec tables; declared states from the GUI-entry oracles),
  `PROT_QC_STATES`. `action-scripts.md` (outside the pin) lists each route's step ids, pinned by
  `test-action-scripts.R`, and four open points for G1(b)/G12.
- `driver.R`: `valGuiHelpers` (extraction from a checkout, `install_enrichment_doubles` as the no-op in
  the extraction parent), `valFixtureCell` (a route on its e2e fixture), `valRunWorkflow` /
  `valWorkflowCell` (the frame of `guiEntryRunLane`, the browser double extended by
  `valInstallDoubleExtension`, rendered `uiOutput`s and modals read after each action, ledger before
  each action, action log for parity, step outcomes A21, forced options and the arm's options scoped by
  `valWithForcedOptions`, the annotation rule `valAnnotationRule`, standard tables then `on_complete`),
  `valReadArmTables` (DA tables read from the route's own written files through `valDaMap`; matrices
  from the state objects by the stage rule), `valReplay`, `valMaskLedger`, `valMaskedDigest`; one
  small fix to group (c)'s `valLedgerSnapshot` (an empty state).
- `arms.R`: `valArmShipped`, `valArmRuvLegacy`, `valArmRuvCalibrated`, `valArmRuvOff`,
  `valArmRuvDesign`, `valArmRuvSupplied`, `valArmInputs` (per-route values read from the panels at
  a60c4e9a, G3 to confirm), `valCellRuv`, `valRunCustomArm`.
- `parity.R`: `valBrowserValues`, `valParityBrowser` (replays the driver's action log through an
  AppDriver), `valCompareLedgers`, `valCompareParityDa`. `integration.R`: `valIntegrationFunctions`,
  `valRunIntegration`.
- `harness/templates/` (the seven files of 7.2), `USAGE.md` (outside the pin).

**Not written (section 7 still open):** `val-run.R`'s cell phases (`lifecycle`, `counted`, `score`,
`determinism`, `correction`, `amended`: steps 5 to 9 of 7.8, wiring `prepare.R`, `valRunCells`,
`valRunWorkflow`, the hook and scoring together; they refuse today); `val-parity.R` and the
`distributable.yml` stage `parity` (7.9: it must be one commit on top of KIT-176's merged workflow,
which this branch, cut from a60c4e9a, does not have; a merge or rebase of `v0.9-staging` is the
coordinator's call); the `run-val-mofa` lane (7.10). Kernel arms beyond the identity gate (a ticket's
kernel function is called by the cell runner, not yet wired).

**Known provisional points in the driver core (to be settled by its first run and by G1(b)/G3):**
the DA column maps (`valDaMap`: metabolomics feature column `metabolite_id` is a reading, not seen); the
display matrix of DIA_LIMPA and transcriptomics cells (read from the state's `args`, which may not hold
it: T-TAB `prot_dia_limpa` will show); `ruv.tsv` reads `args$ruv_factors` and otherwise writes one
`not_recorded` row; the transcriptomics matrices are the route's own transform (TMM, log-CPM, prior
2) recomputed from the QC-complete counts; the DA confirmation modal is assumed to have one primary
button; `MockShinySession` internals are read to find `renderUI` outputs (`private$outs`, as
`guiEntryOutputs` does).

**The next run's exact scope (committed, unrun), for whoever runs the whole suite next:** the whole test
suite at the head (`val-tests.R <out>` with no `--files`), which holds group (a)'s `spec[["target"]]`
fix and the split, group (b)'s `valRunAsChild` append-mode fix, group (c)'s code with the `diffRows`
test fix and the sweep, and group (d)'s code with its tests (incl. the two `test-arm-tables.R`
corrections of a0685450). Expected: groups (a) to (c) passing as stated before; group (d) is new code
never run, so failures there are likely and are first-run findings. T-TAB runs the ten routes end to
end (about 10 to 40 minutes in all; set `timeout_minutes` to 90 and `file_timeout_minutes` accordingly).
Under ruling 1 group (d) has one correction pass of its own.

**After the owner's ruling:** that run; then T2 from `t2/t2-plan.tsv` (three requests A, B, C;
`val-tests.R <out> --coverage-plan docs/qa/workflow-execution/kit/kit-159/t2/t2-plan.tsv --request A
--nsim 2000 --chunks 4`, script limit 280 minutes, combine the three `coverage-table.tsv` files); then
T1 parts with `covr` and `lintr` over the whole package; then R1 (record twin of one small open https
file, library table, phase `pre`; no live UniProt query for the `prot_dia` fixture's synthetic
accessions, ruling 4; the live incomplete-answer refusal waits for KIT-175); then the cell phases of
`val-run.R`; then S0; stop before the release commit (G1(a), G1(r), L1a and later, P1).
`prot_dia_limpa` awaits D19.

**Run-request mechanics:** `.github/run-tests/new-request.sh <name>` in
`/home/doktersmol/Projects/APAF/MultiScholaR-run-tests`; **read the new blank file before writing it**
(a write without the read failed once and the blank file was committed and refused by the plan job:
T0d-blank in `runs.tsv`); a script's `out` comes back in the run's `scripts` artifact (`gh run download
<id> -n scripts`), not `results`; `covr` and `lintr` need `extra_packages`.

## Rulings of 2026-11-23

Recorded as the coordinator gave them (engineering rulings, for the owner):

1. **T1b correction run: GO, once.** Reason: the ticket's one correction pass was written for T1 as one run; because nothing runs on this machine T0 and T1 were split into three module groups (scoring; records and the licence guard; the driver), each shown failing first on a runner. Each group is its own bounded verification item, so the rule is applied per group: ONE correction pass for group (a) (used), ONE for group (b) (this run), ONE for group (c). T1's ceiling of 6 core-hours is unchanged and all three groups and their corrections must fit inside it. A group that still fails after its one correction stops and comes back; nothing is renamed to get another. This is stated now, before group (c) starts, and it applies to instrument defects (a test or script fault) and code defects alike.
2. **T-LOG:** the ticket's own sentinel uses the same `file("/dev/stderr")` call that truncated the private stderr file. Fix the test, and correct the ticket text's sentinel so that it writes the descriptor-2 line first (revision 1.10, closure row); a correction of the instrument, not of a rule.
3. **T2:** splitting each coverage cell's simulations into seeded ranges is approved as implementation, on these conditions: the split is by a declared, fixed seed range per chunk (no seed depends on the chunk count); a test shows the split result is IDENTICAL to the unsplit one on a small cell (same simulated data sets, same counts, same coverage) and fails if a range is dropped or run twice; T2 is sized on the measured factor 2.19 as section 11 says (a third request, +20 to T2's ceiling). No cell may run past the 280-minute script limit by design: ranges sized from the probe with section 11's margin; the probe's per-cell figures recorded. If the projected T2 total passes its ceiling after the third request, stop before sending and report the figure.
4. **R1 and KIT-175:** KIT-175 is blocked on an owner decision. Do everything of R1 that does not need KIT-175's classed UniProt conditions; for the incomplete-answer refusal keep the mocked test T-ANN-1; the live part waits for KIT-175's merge (named open item). Do NOT send a live UniProt query for the `prot_dia` fixture's synthetic accessions (P00001 to P00400); if R1 as written does that, stop that part and report what the ticket says there.
5. The six revision 1.8 additions not yet built belong to group (c) and after, with their tests.

**Open item (ruling 4): R1's live incomplete-answer refusal (G6 in R1) waits for KIT-175's merge.**

## Log

### Step 0 (2026-11-23)

The third re-check's note applied in the main checkout's `tickets/kit/KIT-159.md`: 2.1's heading
"**Decision rule** (the register's, unchanged)" replaced by "**Decision rule** (the round-2 re-check's,
with re-check2 C13's conditions)"; revision 1.7; closure section 22; history line. The rule's text is
unchanged.

### G1D recorded (2026-11-23; not run by this agent)

Run 45121359461, script commit c9c4a412 (branch `run-k159-g1d`), about 1.1 core-hours. Tables and
script copied to `g1d/` (aggregates, column names and digests only, checked before commit); entry in
`gates.md` under G1(r). Observed: the same DA table on 28d5a4be and 12d68ca7 on the runner; head
reproduction holds, base reproduction fails (8a4c2dd5 on the runner against the oracle's 80cfedaa),
so `limpa_da_last_digits` is not admitted on G1D alone. `prot_dia_limpa`: awaiting D19 (owner).

### T0a (group a tests only): run 45122558749 at 525f9733

Observed: `msrval` (no R code) installs; all 53 group (a) tests error with "could not find function"
(`t0/t0a-harness-tests.tsv`, `t0/t0a-harness-failures.txt`). Scripts job 13.2 min, 13.7 job-minutes
with the results job: 0.91 core-hours of T0's 2 expected / 6 ceiling. The script's `out` comes back in
the run's `scripts` artifact (not `results`): `gh run download <id> -n scripts`.

Test corrected after T0a (recorded): T-M3 now passes `n_min = 1`, as 10.1's preamble states for every
hand case (the T0a version omitted it, so the 10-feature floor would have made the value NA).
`test-coverage.R` (pass bands, the regime summary, small regime runs) was written with the coverage
code and has no failing-first run.

### T1a (group a implementation): run 45124309425 at 5ab037ae

Observed: 57 tests, 49 pass (`t1/t1a-harness-tests.tsv`, `t1/t1a-harness-failures.txt`). Seven hand
cases (T-PT, T-JK-1, T-JK-2, T-JK-crit, T-CP, T-ONE, T-REF-FILL) fail with "actual != expected but
don't know how to show the difference": the tests used a relative tolerance (1e-6 or 1e-9) tighter
than the ticket's printed rounding (7 to 10 decimals), so values equal to print precision failed. One
test (small C1 run) recycled a filtered vector against the full level column. Inferred: test defects,
not implementation defects (the waldo message means the difference is below print precision; the
correction run shows the values). `covr` is not in the scripts job's library (it needs
`extra_packages`, as `lintr` does), so the script stopped before covr, lintr and f_sim.

Correction (T1a's one correction pass): the hand-case tests compare within one unit of the ticket's
last printed decimal (`expectClose`, which prints the actual value); the C1 test fixed. covr and lintr
move to the final T1 part over the complete package (A2 is over the whole `msrval`).

T0b (group b tests only, 51ed29f1) was sent before this; its result is recorded when it returns.

### T0b (group b tests only): run 45124916341 at 51ed29f1

Observed: 93 tests; all 36 group (b) tests fail: 23 with "could not find function", the others with
the missing internal constants `HARNESS_PATH` and `EVIDENCE_COLUMNS`; T-DET-1 also had a test defect
(`data.table::fwrite` has no `digits` argument; corrected after T0b, the test now writes with
`fwrite`'s round-trip form). The git-repository fixtures ran (git is in the job). T0 total: 1.46 of 6
core-hours.

### T1a correction (group a): run 45125695623 at 5e1c3202

Observed: all 57 group (a) tests pass, among them every hand case of 10.1, T-VIEW, T-KEY,
T-BLIND-1 and the coverage machinery (`t1/t1a2-harness-tests.tsv`); the 36 group (b) tests fail by
design at this commit. f_sim (G14 (c) on simulations): regime C3, 112 cells at 200 simulations, 94.4 s
in one process against 43.1 s in four concurrent chunks on a 4-core runner: **f_sim = 2.19**
(`t1/t1a2-fsim.tsv`). Cost probe, one cell per family at 20 simulations and 200 target simulations
(`t1/t1a2-cost-probe.tsv`). T1 so far 1.89 of 6 core-hours.

**T2 sizing (inferred from the probe, linear in simulations; not yet a request).** Extrapolated to
2,000 simulations and 20,000 target simulations per cell: about 12 process-hours in all (the ticket
estimated 17.5), dominated by C7 (accuracy about 3.2 h, threshold about 2.1 h, each one cell), C4e at
4,840 features (block count 5: about 2.0 h; 14: not probed, estimated 3 to 5 h) and C4 at 4,840
features (about 0.85 h for g = 4). Two consequences:
1. f_sim = 2.19 is below the ticket's 2.28, so T2 needs a third request (+20 to T2's ceiling, 10.4,
   section 11), which T1 was to show before T2: shown here.
2. A cell is today indivisible (one seeded stream per cell), and at f = 2.19 each of four concurrent
   chunks runs at about 0.55 of solo speed, so a 3-hour cell takes about 5.5 hours of wall and
   overruns the 280-minute script limit. Before T2 the runner must split a cell's simulations (and its
   target simulations) into seeded ranges (one seed per simulation) and combine them; this changes
   `coverage.R`'s seeding and needs its own test (the split result equals the unsplit one). T2 is not
   sent until then.

### T1b (group b implementation): run 45127184692 at c9565ca9 -- STOPPED HERE

Observed: 93 tests, 84 pass (`t1/t1b-harness-tests.tsv`, `t1/t1b-harness-failures.txt`): all 57 of
group (a) and 27 of 36 of group (b) (all of T-IN, T-ANN-1, T-PRE-1, T-CHUNK-1 to 7, T-RUN-1, T-RNG-1,
T-LIC-1 to 11, T-PROV-1, T-DET-2, T-SCORE-1, T-DECL-13, T-DECL-15, T-DECL-17, T-ENV via T-DECL-8 not
reached). Nine fail, each on a test defect read from its message:
- T-DECL-1 to 12 and 14 (seven tests): the test helpers `gitDo`/`gitRepo` pass their arguments to
  `system2` unquoted, so a multi-word commit message (`-m declaration lifecycle`) becomes a pathspec
  and the commit fails; the check then correctly finds no committed declaration
  (`val_declaration_missing`) or the files swept into a later commit (`not_committed_alone`).
  T-DECL-13 passed only because a later one-word commit committed the files.
- T-LOG-1/2: the child script's `writeLines(..., file("/dev/stderr"))` reopens the private stderr file
  in write mode and truncates it, losing the earlier `message()` sentinel. The sentinel did not leak
  (the returned template lines held none); the test's "every sentinel is in the private log" fails.
  Note for the ticket's sentinel (5.2): this call truncates a file-backed stderr, so the sentinel run
  should write the descriptor-2 line first or the check reads it as absent from the log.
- T-DET-1: `data.table::fwrite` writes 15 significant digits, so the 1e-15 change was lost and the
  tables compared `identical`.

Decision under the ticket's attempts rule (section 11: one correction pass for the whole ticket; the
brief: one correction pass per stop rule): T1's correction pass was spent on T1a, so T1b's failures
stop this step and are reported, though all nine are test defects. The three fixes are committed but
not run: `gitDo`/`gitRepo` quote their arguments (`shQuote`); the T-LOG script writes to descriptor 2
first; T-DET-1 writes values with `%.17g`. Needed to continue: the coordinator's word on a T1b
correction run (about 0.6 core-hours; T1 so far 2.5 of 6).

### T1b correction (ruling 1): run 45129042733 at ae464305 -- group (b) STOPPED

Observed: 93 tests, 92 pass (`t1/t1b2-harness-tests.tsv`): the seven T-DECL tests, T-ENV-1 and
T-DET-1 now pass, so those were test-script defects as reported. **T-LOG-1/2 still fails, and this is
a code defect after all:** `valRunAsChild` lets processx open the private stdout and stderr files
without append mode. The child's `writeLines(..., file("/dev/stderr"))` (now first, as revision 1.10
orders) writes through a fresh descriptor at offset 0. The next `message()` on descriptor 2 writes at
that descriptor's own offset 0 and overwrites it (both lines are 18 characters), so the descriptor-2
sentinel is missing from the private log. Nothing reached the template lines or any public file; the
private log is incomplete.

Under ruling 1 group (b) has used its one correction and stops here, to come back. A fix is committed
but not run: `valRunAsChild` starts the child through `sh -c 'exec ... >> out 2>> err'`, so both
private files are in append mode and a reopened `/dev/stderr` cannot overwrite earlier text. It needs
the coordinator's word for a run (group (b) has no correction left).

### Coverage range split (ruling 3): T0s run 45129239148, T1s run 45129513148

- T0s at 64bbbc7c: the three new tests (split identity with data digests, refusal of a dropped or
  repeated range, per-index seeds, the range plan) fail on missing functions; the earlier coverage
  tests pass.
- Implementation at 1a300abe: each simulation, target simulation and C7 bootstrap has a seed fixed by
  the cell's index and its own index only (`valSimSeed`, `valTargetSeed`; ranges 1e9, 1.5e9, 2e9 plus
  index; never the chunk count); `valCoverageRange` runs any range of a cell, `valCombineCoverage`
  joins ranges and refuses a dropped or repeated one; `valCoveragePlan` cuts each cell into a declared
  number of consecutive ranges; `valRunCoverage` runs a plan's request in concurrent chunks (the
  dealing of tasks to chunks changes no result).
- T1s: 93 of 95 pass. **Code defect:** the split-identity test and the small-run test fail because
  `valCombineCoverage` read `spec$target`, which R's `$` partially matches to `target_one` for a cell
  with a simulated target (the function is returned; `notSubsettableError`). Fixed to
  `spec[["target"]]`, committed, not run. Under ruling 1 this code is group (a)'s, whose correction is
  spent, so it waits for the coordinator's word with group (b)'s T-LOG fix. The cost probe of the
  heavy cells came back (`t1/t1s-cost-probe.tsv`).

### T2 plan and sizing (ready, not sent)

`t2/t2-plan.tsv`: 190 tasks over the 169 cells, three requests. Process seconds per cell at 2,000
simulations and 20,000 target simulations, from the two probes (each phase x100; probes at 20 and
200, `t1/t1a2-cost-probe.tsv`, `t1/t1s-cost-probe.tsv`):

| Cell | Estimate (s) | Ranges | Request |
|---|---|---|---|
| C7-159 (bootstrap, accuracy) | 11,556 | 8 | A |
| C7-160 (bootstrap, threshold) | 7,618 | 6 | B |
| C4e-150 (4,840 features, 14 blocks) | 6,459 | 5 | B |
| C4e-148 (4,840 features, 5 blocks) | 7,107 | 5 | C |
| C4-146 (4,840 features, g = 6) | 1,777 | 2 | C |
| C4-144 (4,840 features, g = 4) | 1,268 | 1 | C |
| C4e-149, C4e-147, C4-145, C4-143 (500 features) | 709, 485, 150, 118 | 1 | C |
| C3 (112 cells, each at the dearest probed C3 cell) | 3.5 each | 1 | A |
| C1, C2, C5, C6, C8 | 1 to 2 each | 1 | A |

Per request: A 12,034, B 14,077, C 11,615 process seconds; at f_sim = 2.19 with 10 minutes of set-up,
expected wall 102, 117 and 98 minutes (limit 240 expected, 280 script); the largest task 1,444 s
solo, about 44 minutes at the four-chunk speed of 0.55. Total 10.5 process-hours (the ticket's
estimate was 17.5). Expected 21.1 core-hours; ceiling 60 (40 + 20 for the third request, as section
11 provides); inside it. Every regime row of the coverage table falls within one request (keys carry
the unit count), so the three requests' tables join without recombination. f_sim was measured on C3;
the limma-heavy and bootstrap cells may scale differently, which the margin to 240 minutes covers
about twice over.

T2 waits for the split's fix to pass its run.

### T0c (group c tests only): run 45131922735 at e25055d3

Observed: the five group (c) test files only (`--files`, so no unrun fix of groups (a) and (b) was
exercised): 25 tests, all fail (`t0/t0c-harness-tests.tsv`): missing functions, and three chain tests
first on a test-helper defect (`diffRows` built a zero-row data frame with length-one columns;
corrected after T0c). T0 now 2.84 of 6 core-hours.

### Group (c) code and the defect sweep (committed, unrun)

Written after T0c, not run (the coordinator's instruction of 2026-11-23):
- `driver.R`: `valUiDefaults`, the input state and ledger (`valInputState`, `valApplyRendered`,
  `valApplyServerUpdate`, `valSetInput`, `valScriptSet`, `valLedgerSnapshot`), modals
  (`valApplyModal`, `valModalPress`), `valStepOutcome`, `valRunSteps`, `valWithForcedOptions`,
  `valCheckActionScript`, `valMatrixStages` (revision 1.8 item 4 as a rule), `valExtractHelpers`.
- `gates.R` (a file 7.2 does not list; it holds G1(a)'s chain reader): `valReadChain`,
  `valAdoptReferenceRows`.
- `arms.R`: `valKernelIdentity`, `valCallHook` (named arguments only through named formals, C22),
  `valCopyResultTables` (revision 1.8 item 3).
- `evidence.R`: `valWidenFromRequests` (T-SCORE-2), `valWriteResults`, `valWriteHandoffSha`;
  `licence.R`: the hand-off folder check (`valCheckHandoff`, C21) and the `ticket-coverage.tsv` and
  `ticket-coverage-cells.tsv` columns; `check-evidence.sh`; `val-run.R` with phases `record`, `pre`,
  `check` (the phases that run cells refuse until the driver core exists).

**Not written** (group (c)'s driver core, needed before S0 and the lifecycle): `valRunWorkflow`
(driving `mod_<omic>_server` under `testServer` with the extracted helpers), `valActionScript` and
`action-scripts.md`, `valReadArmTables`, `valReplay`, the RUV arm builders of 7.5, `valRunCustomArm`,
`parity.R` and the `distributable.yml` stage, `integration.R` (Part B), `USAGE.md`, `templates/`,
and T-TAB (which needs the driver core on the fixtures).

Defect sweep (the coordinator's point 3), committed, unrun:
1. Partial matching by `$`: one live defect found (`spec$target`, fixed earlier). Every other `$` on
   a list or record whose element may be absent now uses `[[` with the exact name: `units.R` (eight
   cell-result and condition fields), `driver.R` (`valRunSteps`' step results), `intervals.R`
   (`interval$design`, `interval$transform`, `interval$label`, `accept_flagged$label` and
   `$ledger_row`, and every `member$` field in `valClassifyMember`), `declaration.R`
   (`amend$code`, `$harness_tree`, `$label`; `decl$preparation`, `$hooks`, `$outputs`,
   `$changes_from_lifecycle`, `$arms`; `o$key`, `$labels`, `$max_rows`; provenance `x$stage`),
   `licence.R` (`output$placeholders`, `o$key`, `o$placeholders`, `o$labels`, `o$max_rows`,
   `o$lineage`), `evidence.R` (the provenance fields of `valScoreRuns`), `coverage.R`
   (`regime$kind`, `r$achieved`). Data frame columns and records whose fields are always present
   were left as they are.
2. Files a child process shares for writing: `valRunAsChild` (fixed earlier) and the two chunk
   launchers (`valStartChunks` in `units.R`, `valRunTaskChunks` in `coverage.R`) now start their
   children through one helper, `valRscriptProcess` (`utils.R`), which opens the log in append mode
   (`>>` through `sh`). Note: append mode protects against overwriting, not against truncation; a
   cell that reopens `/dev/stderr` for writing still truncates its chunk log (O_TRUNC), which is why
   the sentinel writes that line first (revision 1.10). `valStatus` already appends; the markers and
   result files are one file per cell.

**The next run's exact scope (committed, unrun):** the whole test suite at the head commit, which
holds: group (a)'s `spec[["target"]]` fix and the split; group (b)'s `valRunAsChild` append-mode
fix; group (c)'s code and the `diffRows` test fix; the sweep. Expected: every test of groups (a),
(b), (c) passing.

### T0d (group d tests only): run 45134811170 at 6d12a614

Observed: the nine group (d) test files only (`--files action-scripts,workflow,arm-builders,
arm-tables,tab,replay,custom-arm,parity,integration`, so no unrun fix of groups (a) to (c) was
exercised): 48 tests, all fail (`t0/t0d-harness-tests.tsv`, `t0/t0d-harness-failures.txt`): 47 on
missing msrval functions, the `action-scripts.md` test on the absent file. No test-script defect showed.
Scripts job 10.7 minutes, 11.2 with the results job: 0.75 core-hours. Before it, the request
`k159-t0d-6d12a614` was pushed blank (my write of the filled file failed because the blank file had not
been read first, and the blank file was committed); the plan job refused it in 6 s (run 45134707220, no
test ran, 0.01 core-hours); the filled request went under a new name. T0 now 3.60 of 6 core-hours.

Test corrections after T0d, found by reading, not shown by the run (a0685450): two `test-arm-tables.R`
values compared whole numbers read back from a TSV (integers) with doubles by `expect_identical`; now
`expect_equal`.

### Group (d) code (committed, unrun): dcc27fb5

As listed under "Resume here". Written after T0d; nothing of it has run.

### Offered regimes from T2 (USAGE.md; provisional ledger row H2) (2026-11-23)

Section 10.2's rule applied as written to `t2/coverage-table.tsv`: a row's recorded status decides
(pass: offered; conservative: offered, so labelled; fail: not offered; flagged: not offered without
`accept_flagged`; not_offered: no interval). No row re-classified. Table in `USAGE.md` (64 rows by
method x class x form x design x transform x n, eight levels each, with the lowest and highest
simulated coverage over the row's conditions and statistics, joined from `t2/t2{a,b,c}-coverage-
detail.tsv`); generator: a 60-line Python join run on this machine (no R), output checked by eye
against the regime table's statuses. Status counts: 337 pass, 35 conservative, 69 fail, 55 flagged,
16 no interval (512). Ledger row `VAL-<merge>` (H2, provisional) in `docs/methods-decision-ledger.md`.

**For the auditor (the rule's wording leaves these unclear; not decided here):**
1. C1 `exact_median` n = 10 at 98.33 % to 99.67 %: `fail` with simulated coverage 0.995 (normal) and
   0.997 (skewed), above the nominal level; achieved level 0.99805. `valCoverageStatus` takes the
   Monte Carlo SE at the achieved level (0.00099; 3 SE = 0.00296), so |0.995 - 0.99805| = 0.00305
   fails by 0.0001; with the SE at the simulated coverage (0.00158, as `coverage-table.tsv`'s `mc_se`
   column prints) it would pass. The ticket says "within 3 Monte Carlo SE of the achieved level" and
   does not say which SE. The interval also under-covers its achieved level, not its nominal one.
2. Rank rows at coverage near or equal to 1: C4e g = 14 `rank` coverage 1.0000 with MC SE 0 and target
   exactly 0 at every level (`conservative` at 95 %, `pass` from 97.5 %, because the upper band
   min(level + 0.04, 1) is 1 from 96 % on); C4 g = 6 and C4e g = 5 `rank` 0.994 to 0.9995 likewise.
   Above 96 % the band cannot mark over-coverage, and a coverage of exactly 1 with a target of 0
   suggests a degenerate statistic (paired pAUC difference to FPR 0.01 identically 0 at 14 blocks?):
   whether these rows license `rank` intervals is for the auditor.
3. C6 `newcombe_paired` 21 pairs at 95 %: `fail` only because the (0.9, 0.9) condition covers 0.9925,
   above the within-run upper band 0.99; the regime table's `coverage` column prints the minimum
   (0.9765), so the table alone does not show why the row fails (the USAGE table prints both ends).
   Offered at 97.5 % and above with coverage to 1.0000 (the band's upper end is 1 there).
4. C8 `clopper_pearson` is `conservative` on every row by construction of `valCoverageStatus` (kind
   `clopper` never returns `pass`); read as "always offered" per 4.4.

### Section 7 code after the merge (2026-11-23): tests first, then code

**Tests only (e7a0c82c):** `test-phases.R` (T-PH-1 to T-PH-7: the cell table, `prepare.R` in an
isolated environment with `universe.tsv` and prepared hashes, the cell runner by arm kind writing
`da.tsv`, a refusal as `val_arm_refused`, the kernel identity gate before kernel cells and the
refusal of every kernel cell when it fails, per-unit metric rows in `metrics.tsv`'s columns and
blinded in the lifecycle, run-level intervals with the refusal rule and `blinded` outcomes, a whole
lifecycle phase writing `runs/<label>/`) and `test-parity-stage.R` (T-PAR-1 the parity request,
T-PAR-2 the bundle's commit against the release commit over `PKG_PATHS`, T-PAR-3 the workflow's
guard on the one admitted upload per platform job, with two mutated copies refused). These test
names are this build's; section 10 lists no test for the cell phases or the stage, so they were
written for the code of 7.8 steps 5 to 9 and 7.9 items 1, 3 and 5. T0e (below) shows them failing.

**Code (2b18b1e3, unrun):**
- `msrval/R/phases.R`: `valTicketEnv`, `valPrepare`, `valPhaseCells`, `valRunCell` (workflow through
  `valRunWorkflow` with the arm builder, forced settings and the post-cell hook; kernel, custom,
  published, baseline), `valKernelGate`, `valUnitMetrics`, `valRunLevel` (jackknife over
  full/deletion units; `paired_t`/`exact_median` over independent/draw units with `valRefusalFill`;
  compound members of each primary family), `valRunPhase` (steps 5 to 9 into `runs/<label>/`).
- `val-run.R`: phases `lifecycle`, `counted`, `correction`, `determinism`, `amended` (declaration
  check, G15 environment check, retrieval, `valRunPhase`) and `score` (`valScoreRuns`, `valRunLevel`).
- `.github/workflows/distributable.yml` (the one workflow change of 7.9, on top of KIT-176's merged
  file): stage `parity` (dispatch only, on `run-val-parity`, input `parity_request`, the request file
  added once and never changed); the steps of stage `full` to the install and the test tools; the
  lifecycle, hosts, lanes, e2e and figure steps skipped in stage parity; a `parity` step (the
  harness of the request's `harness_commit` by `git archive`, `PKG_PATHS` sourced from
  `packaging/package-commit.sh`); after the uninstall, `parity-check` (`check-evidence.sh`); one
  upload `parity-results` per platform job with `if: ... stage == 'parity' &&
  steps.parity-check.outcome == 'success'`, `path: val-out/parity/`, 30 days; the plan job's guard
  admits exactly those two uploads besides the snapshot job's.
- `val-parity.R` (7.9 item 4); `valParityRequest`, `valCheckBundleCommit`, `valCheckParityWorkflow`
  in `parity.R`; `parity-<route>.tsv` and `parity-da.tsv` on both allow-lists (`licence.R`,
  `check-evidence.sh`).
- T2's `coverage-table.tsv` copied unchanged into `msrval/inst/extdata/` (the classifier's
  `valHarnessRegimes()` reads it there; without it no interval can be classified).

**Recorded divergences and readings (for the auditor):** (a) `metrics.tsv` has no comparison
column, so the per-unit rows of an `intersection` or `union` view name their arm set in `view`
(`intersection(a+b)`); (b) a forced setting whose name starts `MultiScholaR.` or `multischolar.` is
set as an option, any other as an input; (c) the post-cell hook receives
`workflow_data$da_analysis_results_list` as the DA results and `paths$integration_dir` as the
Summary folder (read from the product's names, not run); (d) `clopper_pearson` and
`bootstrap_within_run` comparisons over the run's per-unit rows are not wired in `valRunLevel`
(it computes `jackknife`, `paired_t`, `exact_median`), and a within-run interval in the request's
own comparisons (widening, A38) is not computed in a cell phase; (e) `val-parity.R` reads per-route
inputs (`route_inputs`) and L1's digests (`l1`) from the request, fields the request schema of 7.9
does not list by name.

**Not built:** the `run-val-mofa` lane (7.10 sequences it to start with KIT-170; M1 is KIT-170's).

### Ticket revisions 1.8 and 1.9 (2026-11-23)

The coordinator's amendment pass made the ticket 1.8 (register
`AUDIT-KIT-159-plan-amendment-r1.8-20261123.md`); its four minors C21 to C24 applied by me as written
(revision 1.9, closure section 24, history line). Status of the 1.8 additions in the build:

| Addition | Status |
|---|---|
| Phase `pre` result-table types `sd`, `n`, `abs_diff`, `statistic` | built: `valRunPre` checks through `valCheckResultTable`, whose types are 5.2's list (test T-PRE-1 uses `count`; no test yet for the four) |
| `{unit}` over several data sets; a missing expanded file not an error | built in effect: `valExpandOutputPaths` takes the unit ids it is given (the caller passes the union) and `valCheckEvidence` checks only files that exist; no dedicated test yet |
| `ticket-coverage.tsv` fixed columns; `ticket-coverage-cells.tsv` on the allow-list | not yet (allow-list entry and `valTicketRegime`'s column set to change with a test) |
| Hand-off folder in `valCheckEvidence`, `check-evidence.sh`, `val-run.R check` | not yet (tests first) |
| Hook named arguments (C22), standard tables in `val-work/` (C23) | not yet (driver, group c) |
| `prepareInputs()$result_tables` | not yet (`val-run.R`) |
| `matrix_display.tsv` on lipid and metabolite factor-in-model cells | written, unrun (`valReadArmTables`, `valCellRuv`, group d; tested on a hand-made state in `test-arm-tables.R`) |
| `require: "match"` on a signed slope guard | built: `valClassifyMember` takes `require = "match"` on any guard (T-CLS-compound) |
| Phase `score` widening from within-run rows, T-SCORE-2 | not yet (tests first) |
| G1D aggregates under `g1d/`; `prot_dia_limpa` awaiting D19 | done (commit c6f3f0f2) |

No addition changes finished and verified work beyond adding to it.

## Next

0. **Waiting on the coordinator:** a run of the two prepared fixes (group (b): `valRunAsChild` append
   mode for T-LOG; coverage split: `spec[["target"]]`), both code defects found after the groups'
   correction passes. Then T2 from `t2/t2-plan.tsv`, three requests (A, B, C).
1. (Done) The T1b correction run.
2. Before T2: split coverage cells into seeded simulation ranges (one seed per simulation; test: the
   split result equals the unsplit one), then T2 in three requests sized on f_sim = 2.19 (each request
   at most 240 minutes of expected wall).
3. Group (c): driver, arms, helper extraction (T-HELP), the G1(a) chain reader (T-CHAIN), kernel
   identity (T-KER), T-TAB, `valWriteResults`, `val-run.R` (record, pre, check and the other phases),
   `check-evidence.sh`, templates, and the revision 1.8 additions listed above; tests first (T0c),
   then T1c with `covr` and `lintr` (`extra_packages: ["covr", "lintr"]`) over the whole package.
4. R1 (record twin, phase pre, UniProt cache) once `val-run.R record` exists. Its incomplete-answer
   refusal on the real package needs KIT-175's classed conditions, which are not in `v0.9-staging`
   (a60c4e9a has neither `prot_annotation_uniprot_incomplete` nor `_empty`); until KIT-175 merges only
   the mocked unit test T-ANN-1 covers it.
5. S0, then stop before the release commit (G1(a), G1(r), L1a and later, P1).

Former plan for T1a: the group (a) implementation (`utils.R`, `score.R`, `metrics.R`, `intervals.R`, `baseline.R`,
`coverage.R`), tests passing, `covr`, `lintr`, f_sim on regime C3 and a cost probe of every T-COV
family (T2 is sized on these).

## Not checked

- (S0, 2026-11-24) The whole msrval suite: checked later, run 45298108528, all 204 pass. The
  S0 ledgers against G12's census (L1a's job). Determinism of S0's DA tables across jobs (G8).
- (S0) Generated sizes are declared, not tuned: whether a larger input changes any S0 observation.

- Groups (c) and (d) (driver parts and core, arms, parity comparisons, Part B runner), `val-run.R`,
  `check-evidence.sh`: written, never run. T-TAB (the route column maps) not run.
- Not written: `val-run.R`'s cell phases, `val-parity.R`, the `distributable.yml` stage, the MOFA lane.
- `covr` line coverage and `lintr` (A2): not run (deferred to T1c over the whole package).
- T-COV (T2): not run; R1, S0: not run.
- `valBaselineArm` uses `edgeR::calcNormFactors` as the ticket states; edgeR prints that it has been
  renamed `normLibSizes` (a message, no change of result read).
- The 1.8 additions marked "not yet" above.
