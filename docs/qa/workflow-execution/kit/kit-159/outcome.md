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

## Resume here (state 2026-11-24, R1 stopped on a harness defect; read this first)

Ticket r1.13 (VAL-69). Every msrval test has passed (whole suite 45219635014, corrections
45225211841 and 45265533254; lintr 0). Rank re-run R (45227824000) done; USAGE and H2 updated
from it.

**R1 wiring (coordinator's word, 2026-11-24: charged to T1 as unit tests of harness code, VAL-69's
purpose; T0 was exhausted):** tests only 1992d121, run 45264824991 fails first (0.77); code
d6d3c1e1 (`valRecordUniprot`, `val-run.R record --uniprot`, `uniprot.tsv` on both allow-lists), run
45265533254 passes with lintr 0 (0.90).

**R1 stopped (runs 45268141191, 45268145400 at 69b693af; 1.53 of R1's 4):** the open file was
recorded in both twins, equal to r1.12 (1,086 bytes, SHA-256 566d7b3d...) and to each other; phase
pre wrote `ref-qc-shape.tsv` (21 rows, 6 columns) with provenance in both. The UniProt record
stopped before any query result in both: **harness code defect**, `valRecordUniprot` passed the
taxon as text and the product's `getUniprotAnnotationsFull` formats it with `%d`
(`R/func_prot_annotation_uniprot_helpers.R:33`). Fixed (`as.integer(taxon)`), committed, unrun.
R1 is "one run each" with the ticket's single correction pass; whether resending R1 (both
requests) is that correction or an instrument re-run under VAL-54 is the coordinator's call.
Evidence in `r1/r1a/`, `r1/r1b/`.

**Budget:** T0 5.95 of 6; T1 19.09 of 21; T2 26.93 of 54; R1 1.53 of 4.

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

- Groups (c) and (d) (driver parts and core, arms, parity comparisons, Part B runner), `val-run.R`,
  `check-evidence.sh`: written, never run. T-TAB (the route column maps) not run.
- Not written: `val-run.R`'s cell phases, `val-parity.R`, the `distributable.yml` stage, the MOFA lane.
- `covr` line coverage and `lintr` (A2): not run (deferred to T1c over the whole package).
- T-COV (T2): not run; R1, S0: not run.
- `valBaselineArm` uses `edgeR::calcNormFactors` as the ticket states; edgeR prints that it has been
  renamed `normLibSizes` (a message, no change of result read).
- The 1.8 additions marked "not yet" above.
