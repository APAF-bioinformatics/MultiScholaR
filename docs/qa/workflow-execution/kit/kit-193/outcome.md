# KIT-193 outcome (in progress: statements gate passed with corrections S01-S07 applied; waiting for its closure check)

Ticket `tickets/kit/KIT-193.md` r1.2. Branch `kit-193-lean-proofs` from `v0.9-staging` at 51959495. Lean runs on the
laptop by the owner's exception (VAL-95, VAL-103); no R was run anywhere for this stage.

## S0: toolchain and Mathlib (acceptance criterion A1)

Observed on 2026-11-24, laptop, every Lean/lake call pinned to four cores (`taskset -c 0-3`, `LEAN_NUM_THREADS=4`;
Lake 5.0.0 has no `-j`/`--jobs` option, read from `lake --help`; `lean --help` lists `-j, --threads`), one call at a
time. CPU time per call: `cpu-time.tsv` (bash's `time`, user and system time over the call and every child it
waited for; GNU `time` is not installed on the laptop).

- Toolchain `leanprover/lean4:v4.36.0` (`proofs/lean/lean-toolchain`), fetched by elan 4.4.0 on the first `lake`
  call: archive 607,649,028 bytes (`lean-4.36.0-linux.tar.zst`, Content-Length), installed 3.1 GB; 60.1 s.
  `lake --version`: `Lake version 5.0.0-src+611c117 (Lean version 4.36.0)`.
- Mathlib `rev = "v4.36.0"` (`proofs/lean/lakefile.toml`) resolved by `lake update` to
  **fc4f59e577feb1a4956671bc04ac18aaaf8cd842**; `git ls-remote https://github.com/leanprover-community/mathlib4
  refs/tags/v4.36.0` returns the same commit (a lightweight tag). Transitive dependencies are pinned by commit in
  `proofs/lean/lake-manifest.json` (plausible, LeanSearchClient, importGraph, proofwidgets, aesop, Qq, batteries,
  Cli).
- Mathlib's cache: `lake update` ran Mathlib's update hook, which built the `cache` executable (61 small files of
  Mathlib's `Cache/`, not Mathlib itself) and downloaded 9,117 files, 460,999,207 bytes (`~/.cache/mathlib`),
  unpacked to 6,817,999,615 bytes (`.lake/packages/mathlib/.lake/build`). Git checkouts of the nine packages:
  172.6 MB. The whole `lake update` took 98.6 s elapsed. `lake exe cache get` afterwards: "No files to download",
  4.9 s.
- `lake build --no-build Mathlib`: "All targets up-to-date (18249 jobs)": every Mathlib module is served from the
  cache; nothing of Mathlib was compiled on the laptop.
- A trivial root module builds with `lake build --wfail` (third attempt; see below).

Two facts of Mathlib v4.36.0 found at S0 (not in the ticket):
1. `Mathlib.Logic.Basic` no longer exists (first trivial build: "no such file").
2. Mathlib uses Lean's module system. A file that is not a `module` and imports Mathlib gets a Lake warning
   ("designed for use with the module system; consider adding `module`"), which `--wfail` turns into a failure.
   Every project file is therefore a `module` with `public import` and `@[expose] public section` (so the Proof
   files can unfold the Spec definitions). The CI record (4.4) inherits this.

## Statements (section 5.3 item 1)

Spec modules, each built alone with `lake build --wfail MultiScholaRProofs.Spec.<name>` and together through the
root (`lake build --wfail`, 6,136 jobs, all from the cache except the project's 9 modules); no Proof file exists.
`grep -rnE 'sorry|admit|axiom |native_decide|implemented_by|extern'` over `MultiScholaRProofs*`: none.

| File | Statements |
|---|---|
| `Spec/FilterCount.lean` | 8 statements, 6 witnesses |
| `Spec/ReplayRule.lean` | 8 statements, 4 witnesses |
| `Spec/GoClosure.lean` | 12 statements, 4 witnesses |
| `Spec/QvalueBh.lean` | 16 statements, 5 witnesses |
| `Spec/Permutation.lean` | 9 statements, 4 witnesses |
| `Spec/RuvAlgebra.lean` | 12 statements (R-A1 to R-E2, R-K1), 5 witnesses |
| `Spec/RuvDownstream.lean` | 11 statements (R-F1 to R-H1, R-M1), 6 witnesses |
| `Spec/RuvIIIC.lean` | 6 statements (R-J1, R-J2, R-L1, R-L2), 4 witnesses |

Statements gate: `tickets/kit/audits/AUDIT-KIT-193-statements-20261124.md` (Fable). Its corrections S01 to S07
are applied as written (82 statements, 38 witnesses after them). One is a result: the ticket's 4d-ii formula is false
at a = 1 (`kept_family_unit_a_counterexample_stmt`; the formula is now stated for 0 < a < 1).

**Statement not yet frozen:** R-I1, R-J3, R-N1. Each gets its own Spec file and statements review before its proof.

`models.md` holds the model-to-code table and each pair's kind and planted defect. `tests/testdata/kit193/bindings.tsv` lists the bound objects; every digest is `pending` until
its pair first runs on merged code (8.2). `SPEC.sha256` is not written: it follows the gate.

## Freeze

`proofs/lean/SPEC.sha256` written from the Spec files at 4183eb10, equal byte for byte to the values in the
statements register's closure check; committed alone (c783e1ae). From here a Spec file changes only through 5.3 item 4.

## Target 3 (GO closure): proved

Observed on the laptop, 2026-11-24 (four cores, one Lean process at a time). Proof file
`MultiScholaRProofs/Proof/GoClosure.lean`; every statement of `Spec/GoClosure.lean` (12) and every witness (4) has a
theorem whose type is the frozen statement (`Audit.lean`'s `example` lines compile), and each depends on
`propext`, `Classical.choice` and `Quot.sound` only (`#print axioms`; `axioms.tsv`).

- **Laptop CPU time for target 3: 39.8 CPU-seconds** (30.5 s user + 9.3 s system) over 11 Lean calls (the oracle
  module's build and run, six compiles of the proof file, the module build, the root build, the audit), **32.5 s of
  Lean wall time**; labels `t3-*` in `cpu-time.tsv`.
- Agent time: 0.12 h for the proofs (11:41 to 11:48), against 14 expected and 24 ceiling; the oracle table
  (`go-closure.tsv`) 0.07 h before it (`effort.tsv`).
- How the search is proved: a loop invariant (every seen pair is a (term, ancestor) pair; every start pair is seen;
  every seen pair off the frontier has all its successors seen; every seen pair is in terms × nodes) and a potential
  (pairs not yet seen, plus one while the frontier is not empty) that falls by at least one on every pass that does not
  stop; so the pass budget never binds and the end state is exactly the closure.

## Target 2 (replay rule): proved

Proof file `MultiScholaRProofs/Proof/ReplayRule.lean`; every statement of `Spec/ReplayRule.lean` (8) and every
witness (4) proved, each on `propext`, `Classical.choice`, `Quot.sound` or fewer (`axioms.tsv`). The round trip (2c)
is proved from faithfulness (2b) under the option state the replay sets, whose value at a fixed marker is that
marker's recorded rule (a lookup lemma that needs the table's distinct options). The witnesses use the table after
KIT-190 (`missing_value_comparison`, `ruv_control_grouping`, with their tokens). One style linter
(`linter.unusedSectionVars`, an unused instance argument) is off in the file's lemma section; it does not bear on
what is proved.

- Laptop CPU time for target 2: 14.3 CPU-seconds (11.0 user + 3.2 system) over 6 Lean calls, 11.8 s of Lean wall
  time (labels `t2-*`). Agent time 0.05 h (11:48 to 11:51), against 10 expected and 18 ceiling.

## Target 1 (filter count rule): proved

Proof file `MultiScholaRProofs/Proof/FilterCount.lean`; every statement of `Spec/FilterCount.lean` (8) and every
witness (6) proved, each on `propext`, `Classical.choice`, `Quot.sound` (`axioms.tsv`). The exactness within the
lcm bound rests on one integer fact: when `k/n > a/n*`, `k n* - a n` is a positive multiple of `gcd(n, n*)`, hence at
least `gcd(n, n*)`, and `gcd · lcm = n n*` with `lcm < 100000` gives `2 n n* < 200000 (k n* - a n)`, which absorbs the
half-unit rounding slack of `P`. The counterexample "whatever the rounding" forces `P = 312` from `RoundedOK` and
checks the rule at it.

- Laptop CPU time for target 1: 4 calls, 14.2 CPU-seconds (11.9 user + 2.3 system), 10.8 s of Lean wall time (labels `t1-*`). Agent time 0.04 h (11:51 to 11:53), against 10 expected and
  18 ceiling.

## The record job and its lifecycle L0 (acceptance criterion A5)

`.github/workflows/kit-193-lean.yml` (ticket 4.4), branch-filtered to `kit-193-lean-proofs` and
`kit-193-lean-lifecycle`, clean checkout, `actions/checkout` and `leanprover/lean-action` v1.7.1 pinned by commit SHA
(11bd7190…, ba3fdfe5…). Six runs on `kit-193-lean-lifecycle` (`runs.tsv`), 17.9 runner-minutes, about 1.2 core-hours at
the house rule (L0 expected 4.00, ceiling 9.00):

| Push | Run | Observed |
|---|---|---|
| 1, a trivial theorem | 45300586693 (after two instrument re-runs) | every step passed |
| 2, a planted `sorry` | 45301077217 | the build failed under `--wfail` ("declaration uses `sorry`"); the axioms step named `sorryAx`; forbidden words |
| 3, a planted axiom, a Spec byte, an oracle byte | 45302003849 (after one instrument re-run) | axiom-audit (2 declarations use the planted axiom), the Spec hash check (`FilterCount.lean: FAILED`), the axioms step, the oracle comparison (`DIFFERS go-closure.tsv`, `SHA256SUMS` failed), forbidden words |

Three instrument defects, each fixed and re-run under VAL-54 (three re-runs, the maximum for this item):
1. The bundled `leanchecker` of v4.36.0 ends with a segmentation fault (exit 139) when given two or more modules or a
   prefix that covers several; one module per call works (reproduced on the laptop). The job runs `leanchecker -v` once
   per project module (about 7 s each on the runner); the root module, which holds imports only, is not listed.
2. lean-action's nanoda step reads the package name from a `[package]` table that this `lakefile.toml` does not have,
   and would export the whole environment under the package, Mathlib included. nanoda is off: D7 (a) (axiom audit and
   leanchecker) applies; (b) needs the owner (a `lakefile.lean`, and a measured nanoda time over Mathlib).
3. The oracle step lost its failure status in a piped subshell; fixed and checked locally with a stub.

Departures from ticket 4.4, for the auditor: the oracle tables are written by Lean's interpreter (`lake env lean --run
Kit193Oracles.lean <dir>`), not by a compiled `lake exe kit193_oracles`, because Mathlib's cache ships no object files
and a compiled executable would compile Mathlib's C on the laptop and the runner; the interpreter runs the same Spec
definitions. The job's own step runs `leanchecker` (item 1). Module CI, which runs the R impact suite on every push to
any branch, was cancelled on each lifecycle push (no R changes there); the three Module CI runs started by the ticket
branch's first pushes were left to finish. Ticket-branch commits since carry `[skip ci]`; the record L1 runs at the
final commit.

## Paired tests

- Target 3: `tests/testthat/test-kit193-pairs-go.R` (the comparison over the 211 cases of `go-closure.tsv` and the
  planted-defect case: `goAncestorPairs` without the closure). Run P0a (split out of P0 by the coordinator), run
  45302637982 at 8e94b260, R 4.6.1 on ubuntu-22.04: **passed, 638 expectations, 0 failed** (`p0a/tests.tsv`): the
  shipped `goAncestorPairs`, `propagateGoAnnotations` and `buildUniprotGoGeneSets` equal the Lean model on every case,
  and the planted defect is caught. Acceptance criterion A8 for target 3 (base comparison and `R CMD check` stay with
  P0/P3).
- Targets 1 and 2: the oracle tables are committed (`filter-thresholds.tsv`, `replay-*.tsv`); the test files are
  prepared under `pending/` and move to `tests/testthat/` when KIT-187 and KIT-190 merge (they call functions that exist
  only on those branches). Still to write at merge: target 2's builder comparison through `buildWorkflowAnalysisRecord`
  on KIT-187 T8 (iii)'s minimal store.

## D7 ruling (coordinator, 2026-11-24, under VAL-103)

The owner's decision was "a second independent checker if it fits the record job, else axiom audit plus
leanchecker". As built, nanoda does not fit: lean-action's nanoda step reads the package name from a `[package]` table
this `lakefile.toml` does not have, and it would export the whole environment under the package, Mathlib included
(L0 run 45300025617). **nanoda was not run; D7 (a) is final:** the axiom audit plus `leanchecker`, the latter one call
per project module because the bundled `leanchecker` of v4.36.0 crashes (exit 139) on two or more modules (L0 run
45298959712; reproduced on the laptop: one module exit 0, any two exit 139).

## Target 4: `bh_sorted_eq_stmt` is false as frozen (a result; ticket 5.4)

Observed (Lean-checked, `Proof/QvalueBh.lean`, `bh_sorted_eq_counterexample : ¬ bh_sorted_eq_stmt`, on `propext`,
`Classical.choice`, `Quot.sound`): the statement equates R's transcription `bhSorted` with the closed form `bh` for
every `p` over any linearly ordered field, with no hypothesis on the sign of `p`. With negative tied values the later
member of a tie group gets a smaller cumulative minimum (`m/(m - r)` times a negative value grows more negative with
`r`) than the closed form. Instance: `K = ℚ`, `m = 2`, `p = (-1, -1)`: whatever order the sort gives the two tied
indices, the second has `bhSorted = -2` against `bh = -1`.

Inferred (not proved): with `∀ j, 0 ≤ p j` the equality holds (the head of a tie group carries the smallest term).
p-values are never negative, so the route is not affected; the gap is in the statement, written at the gate without
the range 4c states. **The item stops here and goes back to the statements gate (5.3 item 4).** No other target-4
statement uses it (4a to 4c and 4d are stated on the closed form `bh`); work on target 4 stopped at this point by the
coordinator's rule ("stop and report at any statement problem"). Laptop time for this attempt: 11.0 CPU-seconds over 5 Lean calls, 10.2 s of Lean wall time (labels `t4-*`).

## Target 4 (q-value and BH): 4a to 4c and 4d-i proved; freeze amended

- **Freeze amended (5474e512, alone):** statements gate re-review 1 confirmed the counterexample, added `∀ j, 0 ≤ p j`
  to `bh_sorted_eq_stmt`, and froze `bh_sorted_neg_ties_counterexample_stmt` and `bh_sorted_eq_witness_stmt`; one line of
  `SPEC.sha256` replaced. Counts: 83 statements, 39 witnesses. The counterexample statement is proved
  (`bh_sorted_neg_ties_counterexample`). The proof of the corrected `bh_sorted_eq_stmt` waits for the auditor's
  confirmation of 5474e512 (coordinator's rule).
- **Proved** (`Proof/QvalueBh.lean`, each on `propext`, `Classical.choice`, `Quot.sound`): 4a `qflag_contains_bh`,
  `flag_none`; 4b `qflag_extra`; 4c `bh_ge_p`, `bh_mono`, `bh_ties`, `bh_stepup`; 4d-i `grid_mass`,
  `grid_superuniform`, `fdr_n20` (3/40 and 1/20), `fdr_n20_strict` (41/1600), `fdr_n25` (3/50 and 1/25, both
  comparisons); witnesses `qflag_contains_bh_witness`, `qflag_extra_witness`, `bh_unit_witness`, `grid_witness`.
- How 4d-i is proved (as the ticket asks, not by brute force over the outcomes): BH's calls on one and on two values
  are characterised (`callsL_one`, `calls_two`: some value at most `q/2`, or both at most `q`; strictly for the route's
  flag); the sum over the `(N + 1)²` outcomes splits by which p-values are 1 (`fdr2_split`); each grid value's condition
  becomes a whole-number condition on its index (`cond_grid`); the remaining sums are expanded and evaluated by `simp`
  and `norm_num` (each value compiles in under ten seconds).
- Not yet: 4d-ii (`kept_family_formula`, `kept_family_instance`, `kept_family_witness`; second to last by D4) and the
  corrected `bh_sorted_eq_stmt` with its witness.
- Laptop time for target 4 so far: 156.2 CPU-seconds over 24 Lean calls, 103.5 s of Lean wall time (labels `t4-*`, the counterexample's included).

## Target 5 (permutation p-value): proved (5a, 5a', 5b, 5c, 5c')

- **Proved** (`Proof/Permutation.lean`, each on `propext`, `Classical.choice`, `Quot.sound`): 5a' `identity_clause`;
  `perm_tolerance_mono`; 5c `perm_floor`, `perm_floor_reflect`, `bh_floor`, `bh_floor_strict`; 5c'
  `bh_floor_sentence_counterexample`; witnesses `perm_floor_reflect_witness`, `bh_floor_witness`.
- **5a proved** (`perm_finite_core`, witness `perm_witness`; each on `propext`, `Classical.choice`, `Quot.sound`), as
  planned: the stacking map `(τ, π) ↦ (τ, τπ_1, …, τπ_B)` is injective into `Perm^(B+1)` (`stack_injective`; only the
  injection is needed for the bound); the count against index 0 of the stacked tuple (others with a statistic at least
  as large, ties counted) is at most the p-value's count at any `ε ≥ 0` (`cntGe_stack_le`); swapping coordinates 0 and
  `i` shows every index has as many tuples with `1 + count ≤ α(1 + B)` as index 0 (`card_smallAt`); within one tuple at
  most `α(1 + B)` indices have that property (`card_small_count`: the one of least statistic among them has every other
  counted against it); double counting over the `(n!)^(B+1)` tuples gives `(B + 1) #E_0 ≤ (n!)^(B+1) α (B + 1)`.
  Witness: `g = (0, 1)`, `T h = -h 0`, one draw (the swap), `α = 1/2`, `ε = 10⁻⁹`: the draw is not counted, `p = 1/2`.
  Laptop: 7 Lean calls, 29.7 CPU-seconds (24.8 user + 4.9 system), 19.9 s wall (labels `t5a-*`). Agent time 0.10 h
  (12:42 to 12:48), against 20 expected and 30 ceiling (5a and 5a' together).
- **5b proved** (`perm_measure`, witness `perm_measure_witness`; each on `propext`, `Classical.choice`,
  `Quot.sound`): the event is measurable (`measurableSet_permEvent`: for fixed draws `p_ε` is a measurable function of
  the factor, a finite sum of indicators of sets built from `T` and from equalities of coordinates; the draws form a
  countable discrete space); `Measure.prod_apply` writes the probability as the `μ`-integral of the draws' probability;
  invariance of `μ` under each relabelling `τ` (`lintegral_map`) moves the integrand to `f ∘ τ`; summed over the `n!`
  relabellings the integrand is the count of 5a with `g = f` divided by `(n!)^B`, at most `α n!` at every `f`
  (`sum_relabel_le`). Witness: the law with mass 1/2 on each of the constant vectors 0 and 1 in `ℝ^4` (invariant
  because a constant vector is unchanged by every relabelling), `T f = f 0`. Laptop: 5 Lean calls, 40.0 CPU-seconds
  (34.4 user + 5.6 system), 24.5 s wall (labels `t5b-*`; the first call compiled the file before the 5b text was in
  it, an editing slip, not a proof attempt). Agent time 0.07 h (12:48 to 12:52), against 14 expected and 24 ceiling.
- Laptop time for target 5 so far: 31.1 CPU-seconds over 8 Lean calls, 23.5 s of Lean wall time (labels `t5-*`).

## Resume here

Order (coordinator, 2026-11-24): target 5 is complete; next the corrected
`bh_sorted_eq_stmt` with `bh_sorted_eq_witness_stmt` (the auditor confirmed 5474e512); then group R priority 1
(R-A1, R-A2, R-A3, R-C1, R-F1, R-G0, R-G2, R-H1); then stop for the proofs audit (not priority 2, 4d-ii or
`kept_family_witness`). Lean calls go through `taskset -c 0-3` with `LEAN_NUM_THREADS=4`, one at a time, each logged in
`cpu-time.tsv`.

## Not checked

- Targets 4, 5 and group R are not proved yet; nothing was run in R.
- The gate's own reads (RUVIIIC at c534151f, `ruv` 0.9.7.1 as locked, RUVSeq) were not done.
- The CI workflow (4.4) and the lifecycle L0 are not written.
