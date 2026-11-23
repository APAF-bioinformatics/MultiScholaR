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
| `matrix_display.tsv` on lipid and metabolite factor-in-model cells | not yet (`valReadArmTables`, group c) |
| `require: "match"` on a signed slope guard | built: `valClassifyMember` takes `require = "match"` on any guard (T-CLS-compound) |
| Phase `score` widening from within-run rows, T-SCORE-2 | not yet (tests first) |
| G1D aggregates under `g1d/`; `prot_dia_limpa` awaiting D19 | done (commit c6f3f0f2) |

No addition changes finished and verified work beyond adding to it.

## Next

1. On the coordinator's word: a T1b correction run on the commit carrying the three test fixes.
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

- Group (c) (driver, arms, parity, Part B), `val-run.R`, `check-evidence.sh`, templates, `USAGE.md`,
  `action-scripts.md`: not written.
- `covr` line coverage and `lintr` (A2): not run (deferred to T1c over the whole package).
- T-COV (T2): not run; R1, S0: not run.
- `valBaselineArm` uses `edgeR::calcNormFactors` as the ticket states; edgeR prints that it has been
  renamed `normLibSizes` (a message, no change of result read).
- The 1.8 additions marked "not yet" above.
