# KIT-193 outcome (in progress: all 142 frozen statements proved; the last 18 (R-I1, R-J3, R-N1) await their audit)

Ticket `tickets/kit/KIT-193.md` r1.5. Branch `kit-193-lean-proofs` from `v0.9-staging` at 51959495. Lean runs on the
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

**Statement not yet frozen:** R-I1, R-J3, R-N1. Each has its own Spec file; reviewed in
`AUDIT-KIT-193-proofs2-20261124.md` section B, edited as its Q06 to Q08 write and committed for the freeze at 2ab4c189
(Q09); frozen when the auditor confirms that commit by reading (`statements-gate.md`).

`models.md` holds the model-to-code table and each pair's kind and planted defect. `tests/testdata/kit193/bindings.tsv` lists the bound objects; every digest is `pending` until
its pair first runs on merged code (8.2).

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
- Target 1: moved from `pending/` after KIT-187 merged and run (P1a, below): `test-kit193-pairs-filter.R` and
  `test-kit193-pairs-filter-table.R` with `helper-kit193-filter.R`; passed (run 45323167019).
- Target 2: the oracle tables are committed (`replay-*.tsv`); the test file waits under `pending/` and moves to
  `tests/testthat/` when KIT-190 merges (it calls functions that exist only on that branch). Still to write at merge:
  target 2's builder comparison through `buildWorkflowAnalysisRecord` on KIT-187 T8 (iii)'s minimal store.

## D7 ruling (coordinator, 2026-11-24, under VAL-103)

The owner's decision was "a second independent checker if it fits the record job, else axiom audit plus
leanchecker". As built, nanoda does not fit: lean-action's nanoda step reads the package name from a `[package]` table
this `lakefile.toml` does not have, and it would export the whole environment under the package, Mathlib included
(L0 run 45300025617). **nanoda was not run; D7 (a) is final:** the axiom audit plus `leanchecker`, the latter one call
per project module because the bundled `leanchecker` of v4.36.0 crashes (exit 139) on two or more modules (L0 run
45298959712; reproduced on the laptop: one module exit 0, any two exit 139).

## Target 4: `bh_sorted_eq_stmt` is false as frozen (a result; ticket 5.4)

Observed (Lean-checked, `Proof/QvalueBh.lean`, `bh_sorted_eq_counterexample : ¬ bh_sorted_eq_stmt` (superseded by `bh_sorted_neg_ties_counterexample`, below), on `propext`,
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

## Target 4 (q-value and BH): 4a to 4c, 4d-i and the corrected bh_sorted_eq proved; freeze amended

- **Freeze amended (5474e512, with the counter-instance's proof):** statements gate re-review 1 confirmed the counterexample, added `∀ j, 0 ≤ p j`
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
- **Corrected `bh_sorted_eq_stmt` proved** (`bh_sorted_eq`, witness `bh_sorted_eq_witness` with `p = (1, 1, 2)`; each
  on `propext`, `Classical.choice`, `Quot.sound`), after the auditor confirmed the amendment 5474e512 (coordinator,
  2026-11-24). The proof works for whatever order the sort gives tied values: at 0-based position `r` of the decreasing
  order `L` the cumulative minimum is the minimum over `s ≤ r` of `m / (m - s) · p L[s]` (`le_cumMin_iff`); every such
  `L[s]` has `p ≥ p i` and at least `m - s` values at or below it, so its closed-form term is at most the sorted term
  (here `p ≥ 0` is used); conversely, for `j` with `p j ≥ p i` and `c = #{l : p l ≤ p j}`, position `m - c` is at most
  `r` and carries a value at most `p j`, so its sorted term is at most `bhTerm p j`. Laptop: 5 Lean calls, 72.9
  CPU-seconds (67.4 user + 5.5 system), 40.2 s wall (labels `t4-sorted-*`). Agent time 0.07 h (12:52 to 12:56).
- **4d-ii proved; target 4 complete** (`Proof/KeptFamily.lean`: `kept_family_unit_a_counterexample`,
  `kept_family_formula`, `kept_family_instance`, `kept_family_witness`; each on `propext`, `Classical.choice`,
  `Quot.sound`). How: (1) `kept_event_iff`: BH over the terms with `p < 1` calls a term exactly when some `r ≥ 1` has
  at least `r` kept values at most `r q / K` (`bh_stepup_gen` on the kept list; list counts turned into counts over
  `Fin m`). (2) `simes_grid`, Simes' equality on a grid, proved for any finite index type by strong induction on the
  number of values: `n` values uniform on `1..G`, thresholds `λ r` whole wherever at most `G`: the number of outcomes
  with some `r` having at least `r` values at most `λ r` is `min(λ n, G) G^(n-1)`; when `g = λ n < G` the outcomes
  split by the set `S` of values above `g` (`freeze_set` takes out the `(G - g)^|S|` choices the event does not read),
  the rest is the same problem on `T \ S` with `G = g`, and `Σ_k C(n,k) x^k (n - k) y^(n-k-1) = n (x + y)^(n-1)`
  (`binom_deriv`) closes the sum. (3) The false discovery rate splits by the set of kept terms; a set of size `k ≥ 1`
  contributes `a^k (1 - a)^(m - k) min(q/a, 1)` (`GridFine` gives the whole-number thresholds), and the binomial
  theorem sums these to `min(q/a, 1)(1 - (1 - a)^m)`. The instance: every `K ≤ 150` divides `lcmUpTo 150`, so the
  thresholds `(5/27)(i/K)` are the grid points `5 i (L/K) / (27 L)`; `(73/100)^150 ≤ (73/100)^30 < 1/1000` gives
  `> 0.185`. The `a = 1` counterexample by direct evaluation of the 25 outcomes (`m = 2`, `N = 4`, `q = 1/2`).
  Laptop: 12 Lean calls, 52.1 CPU-seconds (labels `t4dii-*`; one, `t4dii-decide-probe`, a probe file in the session
  scratchpad, not in the project, showing `decide` cannot evaluate `callsKept`). Agent time 0.23 h (13:38 to 13:52),
  against 30 expected and 45 ceiling.
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

## Group R priority 1: proved (R-A1, R-A2, R-A3, R-C1, R-F1, R-G0, R-G2, R-H1)

Proof files `MultiScholaRProofs/Proof/RuvAlgebra.lean` and `Proof/RuvDownstream.lean` (both added to the root module).
Every statement of the eight items and their witnesses proved, each on `propext`, `Classical.choice`, `Quot.sound`
(`axioms.tsv`; `Audit.lean`'s `example` lines compile):

| Item | Theorems (statement `<name>_stmt`) | Witness |
|---|---|---|
| R-A1 | `Ruv.A1` | none needed (no hypothesis) |
| R-A2 | `Ruv.A2` | `Ruv.A2_witness` (2 runs, 1 feature, 1 set, `k = 1`, `U = e_1`) |
| R-A3 | `Ruv.A3` | none needed |
| R-C1 | `Ruv.C1` | `Ruv.A2_witness` (shared, as `models.md` lists) |
| R-F1 | `Ruv.F1` | `Ruv.F_witness` (`X1 = e_1`, `X2 = e_2` in `ℝ^3`) |
| R-G0 | `Ruv.G0` | `Ruv.G0_witness` (2 runs in one set, `U = (1, -1)`, `X = M = (1, 1)`, every hypothesis together) |
| R-G2 | `Ruv.G2`, `Ruv.G2_sign` (the two-case instance: `X = e_1`, `C = e_2`, `y = e_3` in `ℝ^4`; ratios 4/3 and 2/3) | `Ruv.F_witness` (shared) |
| R-H1 | `Ruv.H1`, `Ruv.H1_weighted` | `Ruv.H1_witness` (`X = I_2`, first row of `T` = `c = (1, -1)`, weights `(2, 1)`) |

How the two non-routine ones are proved:
- R-C1: with `D = a_c'(a_c a_c')⁻¹ α`, `Wα = Y_c D` and `W_c α = centre(Y_c) D`; column means are linear
  (`colMean_mul`: the mean of `AB` is the mean of `A` times `B`), so restoring adds back exactly the column means of
  `Y_c D`, which is what centring `Y_c` removes. The second centring is the identity because the column means of a
  centred matrix are 0 (also at `n = 0`, where both sides are empty sums).
- R-G0: `U'X = 0` (from `R_M U = U`, `R_M X = 0` and the symmetry of `R_M`), so `U'(Y* - XB) = U'Y* = 0` (R-A2). For the
  rank bound (`rank_bound`): the residuals `R_X V` of any `V` with `U'V = 0` are annihilated by `[X, U]'`; `[X, U]` has
  rank `g + k` because `[(X'X)⁻¹X'; W'R_X]` is a left inverse (`U'W = I`); Sylvester's inequality
  (`Matrix.rank_mul_ge`) then gives `rank(R_X V) + g + k ≤ n`. For the adjusted values, `U'1 = 0` (from `Hone`) gives
  `U'W_c = I` and `U'(Y - W_c α) = 0`, and the same bound.

Two notes for the auditor: R-C1's hypothesis `Hac` and R-G2's `X'X` invertibility are in the frozen statements but not
used by the proofs (both conclusions hold without them; nothing is weakened). In `Ruv.G2_sign` one style linter
(`linter.unusedSimpArgs`, a simp argument not needed in every `fin_cases` branch) is off for that theorem only; it
has no bearing on what is proved.

- Laptop: 11 Lean calls, 67.4 CPU-seconds (57.1 user + 10.3 system), 47.2 s wall (labels `rP1-*`). Agent time 0.13 h
  (12:56 to 13:04) for the eight items together, against 38 expected and 62 ceiling (priority 1); not timed per item.

## Proofs audit (first part) and its minors

`tickets/kit/audits/AUDIT-KIT-193-proofs-20261124.md`: no blocker, no major, eight minors. Applied in one records
commit: P01 (`cpu-time.tsv` row 39's malformed user time marked not recorded; the laptop total below; from row 121 on
the command column shows the command as run: `LEAN_NUM_THREADS=4 taskset -c 0-3 timeout 900
/home/doktersmol/.elan/bin/lake …`), P02 (four stale lines), P04 (`models.md`, notes under the item table), P07 (ticket
r1.4, in the main checkout: interpreter-written oracle tables, `leanchecker` once per module, `nanoda` not run), P08
(the Proof modules' three global instances, `uniformDraws_finite` and the two `Decidable` instances of the grid
section, are now `local instance`; the root rebuilt with `--wfail` and `Audit.lean` re-run, 88 theorems unchanged on
the standard axioms). P03: the coordinator adopted the auditor's two sharper witness statements (a new Spec file, below).
P05 and P06 concern the paired tests.

**Laptop total (implementer, all Lean, from `cpu-time.tsv`):** 122 calls, 786.2 CPU-seconds (616.3 user + 170.0
system, row 39's user time not recorded) at the P01 commit; updated at each later commit below. At b404e62e (R-F2): 134 calls, 835.7 CPU-seconds (657.9 user + 177.8 system). After R-B1 to R-K1: 141 calls, 885.4 CPU-seconds (703.6 user + 181.8 system). After priority 2: 146 calls, 900.4 CPU-seconds (715.5 user + 184.8 system). After R-F3, R-L2, R-M1: 152 calls, 923.3 CPU-seconds (734.8 user + 188.4 system). After 4d-ii: 164 calls, 975.4 CPU-seconds (780.6 user + 194.8 system). After the R-I1, R-J3, R-N1 drafts: 169 calls, 986.0 CPU-seconds (this session: 35 calls, 150.4 CPU-seconds, 88.4 s of Lean wall time). After the freeze edit pass: 174 calls, 1,003.8 CPU-seconds (801.6 user + 202.2 system). After R-I1, R-J3, R-N1: 188 calls, 1,071.4 CPU-seconds (856.6 user + 214.8 system).

## P03: two sharper witnesses (new Spec file; proved)

`Spec/SharperWitnesses.lean` (b8a62ed3, committed alone with its own `SPEC.sha256` line; the two definitions as the
proofs register writes them; the auditor confirms the file at the next statements review). Proved in
`Proof/SharperWitnesses.lean`, each on `propext`, `Classical.choice`, `Quot.sound`:
- `Permutation.perm_measure_moving_witness`: the uniform law on the four unit vectors of `ℝ^4` (invariant because a
  relabelling permutes the unit vectors: `e_i ∘ σ = e_{σ⁻¹ i}`), `T f = f 0`, the atom `e_0` and the swap of samples 0
  and 1 (it moves `e_0` to `e_1` and changes `T` from 1 to 0).
- `Ruv.G0_interior_witness`: three runs in one set, `X = M = (1, 1, 1)'`, `U = (1, -1, 0)'`, two features, the first
  the only control; `n = 3 > g + k = 2`; the corrected second feature is `(0, 0, 1)` and its residual on `X` is
  `(-1/3, -1/3, 2/3) ≠ 0`.
One style linter (`linter.unusedSimpArgs`) is off in that file (simp arguments not needed in every `fin_cases`
branch); no bearing on what is proved. Laptop: 8 Lean calls, 29.1 CPU-seconds (labels `p03-*`). Agent time 0.13 h.

## Group R priority 2

- **R-F2 proved** (`Ruv.F2`, witness `Ruv.F_witness`; standard axioms). The Spec's one hypothesis, `[X1, X2]` of full
  column rank, gives both invertibility lemmas, proved, not assumed: an invertible Gram matrix makes `mulVec` injective
  and conversely (`mulVec_inj_of_gram`, `gram_unit`), so `X2` and `R_2 X1` (`R_2 = residop X2`) are injective, the
  latter because `R_2 X1 v = 0` puts `X1 v` in `col(X2)` and `[X1, X2]` is injective. Then `y = X1 β1 + X2 β2 + r` with
  `X1'r = X2'r = 0`, so `R_2 y - R_2 X1 β1 = r` and `(R_2 X1)'r = X1'r = 0`: `β1` solves the normal equations of the
  second fit (`betaHat_of_normal`), and its residual is `r`. Agent time 0.03 h (13:20 to 13:22; labels `rF2-*`).
- **R-B1, R-B2, R-D1, R-D2, R-E1, R-E2, R-K1 proved** (`Proof/RuvAlgebra2.lean`: `Ruv.B1`, `B2`, `D1`, `D2`, `E1`, `E2`,
  `K1`, `K1_gls`; witnesses `B1_witness`, `D_witness`, `K1_witness`, `K1_gls_witness`; each on `propext`,
  `Classical.choice`, `Quot.sound`). How: a shift `1c'` is killed by `R_M` (from `Hone`) and by `U'` (`U'1 = 0` from
  `HU` and `Hone`), so `α` and `a_c` are unchanged and `W` gains `1 c_ctl' A`; the corrected values gain `1 d'`
  (`shifted_Ystar`), which restoring removes up to `1c'` and which every contrast annihilates; E1 is the constant case.
  D1/D2: `α ↦ tQ'α`, `(t²Q'GQ)⁻¹ = t⁻²Q⁻¹G⁻¹Q'⁻¹`, so `W ↦ WQ'⁻¹` and `Wα` unchanged (D1) or scaled by `t` (D2).
  K1 (ii): `Y* = (I - SUG⁻¹U')Y` with `U'SU = G = a_c a_c'`; `K = I - SUG⁻¹U' - M(M'S⁻¹M)⁻¹M'S⁻¹` kills `M` and
  `SU` (`U'M = 0` from `UU' = R_M` and `R_M M = 0`), and `I = M A + S U B` with `B = G⁻¹U'`,
  `A = (M'M)⁻¹M'(I - SUB)`, so `K = 0`. Instances: R-E2 and the B and D witnesses on four runs in one set, `k = 1`,
  `U = (1, -1, 0, 0)'`, one control (offset amounts 0 and 1 on features 0 and 1); the K1 witnesses on four runs in
  one set, `k = 3`, `U` the three non-constant columns of the 4 × 4 Hadamard matrix over 2 (`UU' = I - 11'/4`), all
  four features controls with `Y = I` (so `S = I`). R-K1 (ii)'s hypotheses `IsUnit (M'S⁻¹M).det` and `Hac` are used
  as stated (not derived). Laptop: 7 Lean calls, 49.7 CPU-seconds (45.7 user + 4.0 system), 23.3 s wall (labels
  `rP2-alg-*`). Agent time 0.13 h (13:24 to 13:32), against R-B1 4/6, R-B2 2/4, R-D1 3/5, R-D2 5/9, R-E1 4/7,
  R-E2 2/3, R-K1 7/12 expected/ceiling.
- **R-G1, R-J1, R-J2, R-L1 proved; group R priority 2 complete** (`Proof/RuvDownstream2.lean`: `Ruv.G1`,
  `G1_witness`; `Proof/RuvIIIC.lean`: `J1`, `J1_witness`, `J2`, `J2_contribution`, `J2_witness`, `L1`,
  `L1_invariance`, `L1_invariance_witness`; each on `propext`, `Classical.choice`, `Quot.sound`). R-G1 from R-F2 with
  `X2 = C` (`R_C X = X` because `C'X = 0`; `y* = R_C y`), the Gram matrix of `[X, C]` block diagonal, and the real
  identity `N/√(ρv/a) = √(a/b) · N/√(ρv/b)` for `a, b > 0`, which holds also where a square-root argument is not
  positive (both sides Lean's `x/0 = 0`), as the statement's docstring says. R-J1: on a complete feature `R_j` is every
  run, the mask is `I`, and `W_j`, `a_j` are those of `complete`. R-J2 by `Σ(x_i - c)² = Σx_i² - 2cΣx_i + mc²` at the
  two means; its link through `(1'1)⁻¹ = 1/(m+1)`. R-L1 by `λ'(Y_C w) = (Y_C'λ)'w`; the invariance because `U'1 = 0`
  kills constants in `a` and in `U'Y_C`, and `Σλ = 0` kills them in `d_C`. Laptop: 5 Lean calls, 15.0 CPU-seconds
  (labels `rP2-G1-*`, `rP2-J-*`, `rP2-GJL-*`). Agent time 0.05 h (13:32 to 13:35), against R-G1 4/6, R-J1 5/8,
  R-J2 3/5, R-L1 5/8.

## Group R priorities 3 and 4 (R-I1, R-J3, R-N1 wait for their statements review)

- **R-F3, R-L2, R-M1 proved** (`Proof/RuvDownstream2.lean`: `Ruv.F3_span`, `F3_fwl`, `F3_fwl_witness`, `M1`,
  `M1_witness`; `Proof/RuvIIIC.lean`: `L2`, `L2_witness`; each on `propext`, `Classical.choice`, `Quot.sound`).
  R-F3: R-F1's shift and R-F2's argument with the weighted normal equations `X'Ω(y - Xβ̂) = 0`; positive
  definiteness of `Ω` turns injectivity of `X2` and of `M_2 X1` into invertible weighted Gram matrices (Mathlib's
  `PosDef.conjTranspose_mul_mul_same`), and `M_2'Ωr = Ωr` because `X2'Ωr = 0` (no symmetry of `Ω` is used). R-L2:
  the quadratic bound makes `G = a_C a_C'` injective, hence invertible; with `z = G⁻¹a`, `‖A a‖² = z'a`,
  `κ z'z ≤ z'a` and Cauchy-Schwarz give `z'a ≤ a'a/κ`, then Cauchy-Schwarz once more. R-M1: `[X, WQ]` and `[X, W]`
  reach the same predictors through `β_W ↦ Q⁻¹β_W` and back; `Matrix.rank` is the dimension of that range. Witnesses:
  `Ω = diag(2, 1, 1)` with `X1 = e_1`, `X2 = e_2`; `a_C = [I_2, 0]`, `κ = 1`; `Q = (2)`. Laptop: 6 Lean calls,
  22.9 CPU-seconds (labels `rP3-*`, `rP34-*`). Agent time 0.05 h (13:35 to 13:38), against R-F3 6/10, R-L2 8/14,
  R-M1 4/6.

- **P1a, run 45316817535 (observed):** the test file reached its 15-minute limit (`timeout`) with no failed
  expectation: 183,822 expectations passed before the limit (gap 1 and 6.1 (i) complete, 6.1 (ii) about 60 percent
  through by the count of its expectations); no part was judged as a whole. The digest script ran (`p1a/bindings-t1.tsv`,
  `p1a/versions.txt`: R 4.6.1, MultiScholaR 0.9.0, digest 0.6.39); `bindings.tsv` rows T1 filled from it. Jobs: tests
  22.7 minutes, scripts 9.0; about 3.45 core-hours. **Instrument defect (VAL-54, re-run 1 of 3):** the failure is in the
  test file's run time, not in a comparison. Fix, with no expected value, comparison or limit changed: the file split
  in two (`test-kit193-pairs-filter.R`: gap 1, 6.1 (i), (iii), (iv), the planted defect; `test-kit193-pairs-filter-table.R`:
  6.1 (ii), (v), gaps 2 to 5), each within the ticket's 15 minutes; 6.1 (ii) and (v) checked per (n*, m) over all its
  rows at once (the same comparisons, 5,490 expectations in place of 221,430: three for each of 1,830 groups); the shared helpers in
  `tests/testthat/helper-kit193-filter.R`.

- **P1a re-run, run 45323167019 at c8123f09 (observed): passed.** `test-kit193-pairs-filter.R` 50,942 expectations
  in 452.5 s, `test-kit193-pairs-filter-table.R` 11,533 in 667.8 s, 0 failed, 0 errors (`p1a-2/`). Every shipped cutoff
  satisfies `RoundedOK` (groups 1 to 316; G up to 40); the decision, the floor form and `allowedByRoundedPercent` equal
  the oracle on all 73,810 rows; two-group designs store the larger group's P; the full helper keeps exactly the
  oracle's rows at the group and the protein level; the counterexample rows behave as the model says; the planted
  defect is caught in `exceedsRoundedPercent` and through the helper. Target 1's pair holds (A8 for target 1). P1a in
  all: about 6.0 core-hours (3.45 + 2.58) against P1's 4.67 ceiling for one job (the ticket's R total is 63.34).

## Notes for the auditor (second implementer session)

- Hypotheses in frozen statements that the proofs do not use (nothing is weakened; each conclusion holds without
  them): `Hset` and `Hac` in R-B1, R-B2, R-E1 (the shift argument needs only `HU` and `Hone`); `Hac` in R-D1 (only
  `Q` invertible is used); `Hset` in R-K1 (i); `0 < m` in `kept_family_formula` (at `m = 0` both sides are 0).
- R-K1 (ii) uses all five of its hypotheses as stated; `IsUnit (M'S⁻¹M).det` and `Hac` are not derived (the Spec
  allows either).
- `Proof/KeptFamily.lean` proves Simes' equality on a grid for any finite index type (`simes_grid`), more general than
  4d-ii needs; it is a lemma, not a frozen statement.
- One unrelated reading, left alone: `getNegCtrlProtAnovaHelper` takes `sort(…)[1:list_size]`
  (`R/func_prot_norm_optimization_helpers.R:582`); at `list_size = 0` R's `1:0` is `c(1, 0)`, which selects one
  control. R-N1's draft models it (`rHead`); whether it is reachable is for the owner (not checked).

## Target 1's paired test (run P1a: P1 split, target 1 only; target 2's pair waits for KIT-190)

- **Ticket 8.2's diff before the first run:** `git diff e6c42d1d d84b6dfd -- R/func_prot_qc_filtering_helpers.R
  R/func_general_parameter_helpers.R R/utils_workflow_analysis_record.R` is empty: the lines `Spec/FilterCount.lean`
  cites (read at e6c42d1d on `kit-187-filter-by-count`) are unchanged on `v0.9-staging` at d84b6dfd (KIT-187 merged at
  6416dbd6). Nothing goes back to the auditor.
- **Merge:** `v0.9-staging` (d84b6dfd) merged into `kit-193-lean-proofs` (147d498d, a true merge, clean), as ticket
  11.1's P1 row says ("the branch merged with staging first").
- **The test** `tests/testthat/test-kit193-pairs-filter.R` (moved from `pending/`): 6.1 (i) to (v) and the five gaps
  of AUDIT-KIT-187-impl I10: (1) the marker's option `MultiScholaR.missing_filter_rule` unset by
  `withr::local_options` in every test and asserted unset; (2) `removeRowsWithMissingValuesPercentHelper` driven on 57
  oracle rows (n* in {3, 7, 20, 33, 60}, m in {1, ⌈n*/2⌉, n*}, n in {1, ⌈n*/3⌉, n* - 1, n*}; two groups, the larger
  sets the cutoff, min_groups 2), each with rows k = 0..n missing in the group under test, kept set equal to
  `k0..kmax`; (3) two-group designs a < b ≤ 60 (a in {1, b %/% 2, b - 1}, every m): the stored cutoff's P equals the
  model's P of the larger group; (4) the grouping column `group` in `updateMissingValueParameters` and in the helper;
  (5) the group level through the helper: G groups of one run (G = 1..20 and 40, every g), rows with f groups failing,
  kept set equal to `f0..kmax(G, g, G)`. Planted defect (ticket 5.5): `roundedPercentNumerator` without the half step,
  caught both in `exceedsRoundedPercent` and through the helper (`expect_failure`), after both comparisons are shown to
  pass on the shipped code. P05 of the proofs register concerns targets 5 and R-G0's pairs only; none applies here.
- **Bindings (8.2):** a script, `scripts/bindings-t1.R`, prints `digest::digest(deparse(body(fn)), "sha256")` of the
  four T1 functions on the installed package; its values fill `bindings.tsv` rows T1 after the run.
- **Declared before sending:** purpose: tie target 1's proved model to the shipped code (A8 for target 1). Commit: the
  commit that adds the test (on `kit-193-lean-proofs`, after the merge). Ceiling: one test job of 60 minutes (P1's),
  the file limited to 15 minutes (ticket 6.1); one scripts job, 30 minutes; at most 6.0 core-hours, expected about 3.5.
  Decides: every comparison passing and the planted case caught → target 1's pair holds; any comparison failing on the
  shipped code → a result: target 1's item stops, the owner and KIT-187's coordinator are told (5.5); the planted case
  not caught, a timeout or a crash with no comparison judged → an instrument defect (VAL-54, at most three re-runs).

## Second proofs audit and the statements review of R-I1, R-J3, R-N1

`tickets/kit/audits/AUDIT-KIT-193-proofs2-20261124.md`: the 124 theorems sound (A); the three drafts true and none
to wait, each frozen after one edit pass (B); target 1's pair holds (C). Applied:
- **Edit pass and freeze (Q06 to Q09), commit 2ab4c189**, holding only the three Spec files and their three replaced
  `SPEC.sha256` lines (the nine frozen files and lines unchanged; `sha256sum -c` twelve OK). Q06: R-I1's header (runs
  removed after RUV; KIT-184 r1.4), `I1_UtWc`'s `Hset` sentence, `I1_witness` away from the corners, `I1_rank` added;
  the two optional statements `I1_ctl_rank` and `I1_centred_diff_rows` included (the register's texts are complete),
  so the header's last clause reads as the register gives it for (b'). Q07: `J3_counts`' docstring (the missing-value
  part of the clean-up's count only), the two header sentences, `J3_witness` with `0 < k` and a general-branch feature.
  Q08: the line-580 quotation, which rule is modelled (`q_table_order`; KIT-190's other rule not modelled), the
  row-name and missing-q sentence, the two-part `N1_rules_differ`; the optional `N1_ties` included. The header of each
  file lost its "Not frozen" paragraph and the title its "draft for its statements review". One departure from the
  register's text, layout only: in R-I1's header the paragraph after the inserted sentence is re-wrapped to the file's
  line width. Every file compiled alone (`lake env lean`) and the root built with `lake build --wfail` (7,374 jobs,
  no warning). The auditor confirms 2ab4c189 by reading; no proof of the three is written before that.
- **Records (Q02, Q03, Q05, Q10 (a)):** `models.md` (the unused hypotheses and the corner witnesses; what 4d-ii may be
  cited for; the stale run sentence; T1's two pair files and fifth bound object; section 3 for the frozen texts and
  R-N1's pair; section 4's digest line); `bindings.tsv` (T1 rows name both pair files; a fifth T1 row,
  `removeRowsWithMissingValuesPercentHelper`, digest `pending`); `statements-gate.md` written; the stale lines here.
- **Needs a runner job, not run (Q03 item 7, Q10):** the digests of target 3's five bound objects and of T1's helper,
  and **the shipped P at (320, 137)** (ticket 6.1 (iii): the test asserts only that it is 57187 or 57188). The script
  `scripts/bindings-t3.R` writes both (`bindings-t3.tsv`, `p-320-137.tsv`); it is to go with the next R run. Until
  then those rows stay `pending` and the P is not recorded.
- Laptop: 5 Lean calls, 17.8 CPU-seconds (13.2 user + 4.6 system), 13.4 s wall (labels `freeze-*`).

## R-I1, R-J3, R-N1 proved (statements frozen at 2ab4c189, confirmed by the auditor)

- **R-I1 proved** (`Proof/RemovedFactors.lean`: `Ruv.I1_UtWc`, `I1_centred_diff`, `I1_centred_diff_rows`, `I1_span`,
  `I1_rank`, `I1_ctl_rank`, `I1_average`; witnesses `I1_witness`, `I1_average_witness`; each on `propext`,
  `Classical.choice`, `Quot.sound`). How: centring commutes with multiplication on the right (`centre_mul`, from
  `colMean_mul`), so `W_c = centre W`, `centre(Wα) = W_c α` and `centre(W_c α) = W_c α`; on a selection of runs the
  rows of `W_c` are those of `W` less one constant row, which centring removes (`centre_sub_rows`). `U'1 = 0`
  (`Ut_one`) gives `U' centre(Y_c) = a_c`, hence `U'W_c = I` and rank `k` (Sylvester's bound through
  `rank_mul_le_right`); the same left inverse gives `rank(W_c α_F) = rank α_F`. The spans: `mulVecLin` of a product is
  a composition, and `rank α_F = k` makes `α_F` onto. (c''): `a_c = α_F S` for a selection matrix `S` when
  `ctl ⊆ F`, and `rank(a_c a_c') = k`. (d): `W = W_c + 1 m'`, so `A` sends a removed column to
  `(A W_c) α_j + t 1` (`fromCols_mulVec_sumElim`). Witness: `witnessE2` with `F = {0}` (`W = e_1`,
  `W_c = e_1 - 1/4`; the rank clause is `I1_ctl_rank` at `F = ctl`); the averaging map `(1/2, 1/2)`. Laptop: 5 Lean
  calls, 32.9 CPU-seconds (labels `rI1-*`). Agent time 0.10 h, against 10 expected and 16 ceiling.
- **R-J3 proved** (`Proof/RuvIIICPattern.lean`: `Ruv.J3_pattern`, `J3_counts`; witness `J3_witness`; each on
  `propext`, `Classical.choice`, `Quot.sound`). Both by unfolding `out` per run (a kept feature is a corrected one, and
  then a value is missing exactly where the run is not observed; `Finset.filter_congr`). Witness `witnessJ3`: four
  runs in one set, `k = 1`, feature 0 (the only control) observed on runs 0 to 2 (residual dimension 2 > k, the
  general branch, the eigen step succeeding), feature 1 observed on run 0 only (residual dimension 0: not corrected);
  the residual dimension is the code's, observed runs minus one set. Laptop: 3 Lean calls, 10.7 CPU-seconds (labels
  `rJ3-*`). Agent time 0.03 h, against 2 expected and 3 ceiling.
- **R-N1 proved** (`Proof/RuvControlRules.lean`: `RuvControls.N1_stable`, `N1_ties`, `N1_first_pass`, `N1_trx`,
  `N1_rules_differ`; witness `N1_stable_witness`; each on `propext`, `Classical.choice`, `Quot.sound`). How: (a) with
  every `q` equal the comparator always holds, so the kept list is already sorted and `mergeSort` returns it (core's
  `mergeSort_of_pairwise`). (a') `mergeSort` output is sorted by "larger q first" (`pairwise_mergeSort`); for two
  tied values in the wrong row order, the filtered list (strictly increasing rows) has them in row order, stability
  (`pair_sublist_mergeSort`) keeps that order in the output, and a list without repeats cannot hold both orders.
  (b), (c) from one counting lemma (`card_rank_gt`): under a strict total order on a set of `N` elements the rank
  `1 + #{j before i}` is a bijection onto `1..N` (injective by totality and `rank_lt_rank`; into `1..N`; equal
  cardinalities), so exactly `N - m` ranks exceed `m`; "smaller p, then earlier row" on the present p-values (`F/4 <
  r` in ℚ is `⌊F/4⌋ < r`), and R's `order()` (missing last) on all genes; the second parts from `rank_lt_rank` and
  totality. (d) part 1 from the counts (4 against 3 controls at five genes); part 2 by evaluation (gene 0 has no
  p-value, so it is not a first-pass control; one present gene precedes it in `order()`, so its position exceeds
  `⌈4/4⌉ = 1`). One limit of Lean met on the way: `decide` cannot evaluate statements that unfold `before` (its match
  compiles to a private splitter the kernel cannot reduce), so (d) is proved by the counts and by hand, not by
  `decide`. Laptop: 6 Lean calls, 24.0 CPU-seconds (labels `rN1-*`). Agent time 0.07 h, against 2 expected and 4
  ceiling.
- **All 142 statement constants proved.** `Audit.lean` has 142 `example` lines (each theorem's type is its frozen
  statement) and 142 `#print axioms` lines, every one on the standard axioms or fewer; `lake build --wfail` clean
  (7,380 jobs); `sha256sum -c SPEC.sha256` twelve OK, no Spec file changed. Stage total: 14 Lean calls, 67.5
  CPU-seconds; agent time 0.20 h against 14 expected and 23 ceiling for the three items.

## P2 and P0b: target 5's pair, target 4's library parts, group R's pairs on merged code (run 45335802998)

Fourth implementer session, 2026-11-24. Branch merged with `v0.9-staging` at ea93d948 first (84e55b0d, a true merge,
no conflict: KIT-175, KIT-177, KIT-189, KIT-191, KIT-199 and the ledger to VAL-137).

**Ticket 8.2's diff before the first run (observed, by `git show` and `awk` over the function bodies).**
`integrationEnrichmentSetT` and `INTEGRATION_ENRICHMENT_PERMUTATION_TIE` (1e-9) are as at 2f54469e, the commit the
Permutation Spec cites; `integrationEnrichmentPermutation` differs only in the values of the time-limit refusal's
message (outside the model). The floor: at 2f54469e the module computed `ceiling(m p_min / 0.05 - 1e-9)`; merged
KIT-191 computes it in whole numbers, `integrationEnrichmentFloorMinimum()`, `(20 m (1 + d) + B) %/% (1 + B)`, which
is `bh_floor_stmt`'s `ceil(m p_min / q)` exactly; the model covers it (the pair checks it on 1.44 million cases).
Group R: `restoreRuvFeatureLevels`, `computeRuvCovariates`, `computeProtRuvModelFactors`, `protRuvFirstPassControls`,
`computeOmicsRuvAssayFactors`, `omicsRuvFirstPassControls`, `ruvExactContrastsFit`, `trxRuvDesignControls` and the
`ruvIII_C_Varying` method for `ProteinQuantitativeData` are byte-identical to 147adde8 (and `trxRuvDesignControls` to
19997323). Nothing goes back to the auditor.

**Oracle tables (Lean, f3170c03).** `Oracle/Permutation.lean` (`perm-floor.tsv`: 3 fixtures x 13,824 cases, `d`, `r`,
`(1 + d)/3`), `Oracle/QvalueBh.lean` (`bh-counterexample.tsv`: 441 + 676 outcomes, the Spec's `callsKept` and
`callsWhole` at "at most" and "below" q), `Oracle/RuvControlRules.lean` (`first-pass.tsv`: 364 cases, the Spec's
`firstPassControls` and `trxFirstPassControls`). Regenerated with `lake env lean --run Kit193Oracles.lean`; the six
earlier tables came out byte-identical; `SHA256SUMS` rewritten; `tests/testdata/kit193/` is 3.0 MB. Checked on the
laptop from the tables (Python, exact fractions): FDRs 3/40, 1/20, 41/1600 (N = 20) and 3/50, 1/25 for both
comparisons (N = 25), as `fdr_n20_stmt`, `fdr_n20_strict_stmt`, `fdr_n25_stmt` state. No Spec, Proof or Audit file
changed. **Owed (not mine to write under this session's brief):** ticket 5.5 asks of a separate executable oracle a
theorem in the Proof files, listed in `Audit.lean`, that it equals the Spec definition. The Permutation Spec's
`unchangedDraws` and `reflectedDraws` use classical decidability (stated for any type of values) and cannot be
evaluated, so `Oracle/Permutation.lean` counts the same filters over lists of rationals with their decidable equality.
The bridge lemma (a filter does not depend on its decidability instance, `Finset.filter_congr_decidable`, and a
function on `Fin 4` is its list of values) is not yet in `Proof/`. The other two oracles evaluate Spec definitions
directly. Laptop: 3 Lean calls, 15.5 CPU-seconds (labels `p2-oracle-*`).

**Pairs (tests only, d71dbd08), one shard each, file limit 15 minutes; all passed, 0 failed.**
- `test-kit193-pairs-permutation.R` (target 5; 30 expectations, 55.7 s): all 576 draw pairs reached through the seed;
  fixture A's 24 values of |t| more than 1e-8 apart and fixture C's 12 (condition met); on A and B `identity_draws`
  = `d` and `min_p` = `(1 + d)/3` on all 27,648 cases; every p at least `(1 + d)/3`; on C `identity_draws` = `d` and
  every p at least `(1 + d + r)/3` (the reflection ties counted); `#{p <= k/3} <= 4608 k` for k = 1, 2, 3 on A and C;
  the floor 60 at B = 1,000, d = 10, m = 273 (the sentence prints 60) and the shipped floor equal to the exact ceiling
  on m 1 to 6,000, d 0 to 119, B in {1,000, 10,000}; through `integrationEnrichmentKernel` (KEGG, mapped list
  mocked), `tests` = the rows `p.adjust` adjusted over (6 of 7 pathways; the one below the size limit is outside the
  family) and the sentence built with that m; 5c' by `stats::p.adjust` (4 called at 50/1001, floor 2, none at or
  below 1/1000). Calls go through `integrationEnrichmentPathwayTest()`, the entry point the KEGG, Reactome and kernel
  routes call; the exported `integrationEnrichFactor()` needs a trained project and was not used. The tolerance at a
  large |t| (proofs audit P05 item 1) is a note, not tested.
- `test-kit193-pairs-qvalue.R` (target 4 library parts; 161 expectations, 7.7 s): (i) on 20 seeded vectors (length 1
  to 200, ties, 0 and 1) and pi0 in {1, 0.5, 0.13}, `qvalue::qvalue` equals `pi0 * p.adjust(p, "BH")` to relative
  1e-14 and nests the BH flag (no exception); the same with qvalue's own pi0 on 5 vectors of 200; (iii) `p.adjust` on
  the 441 and 676 outcomes calls exactly what the model calls, every outcome and both comparisons; FDRs 0.075,
  41/1600, 0.05 (N = 20) and 0.06, 0.04 for both comparisons (N = 25), to 1e-12. Not here: (ii) and the
  `processEnrichments` flag, which KIT-203 changes (binding row T4 `processEnrichments` stays `pending`; ticket 8.1).
- `test-kit193-pairs-ruv.R` (group R; 181 expectations, 11.2 s): R-A1, R-A2 (`ruv::RUVIII`), R-A3, R-B1, R-B2, R-E1
  (with E2's differing amounts on the shipped output), R-C1, R-D1, R-D2, R-F1, R-F2, R-F3 (weights), R-G0, R-G1 (an
  illustration), R-G2, R-H1 (and weighted), R-J1, R-J2 (upstream's zero-filled equality on every feature, to 1e-8),
  R-J3's pattern, R-K1 (k = 9, 14 controls: equal to `ruv::RUVIII` at k = 9 and at k = NULL, within-set range <=
  1e-8), R-L1, R-L2, R-M1 (edgeR fitted values and F to relative 1e-8), R-N1's first-pass rules (all 364 cases of
  `first-pass.tsv`, all three functions). R-G0's "equal n - g - k" is a property of the fixture (general branch,
  default formula); the theorem gives "at most" (P05 item 4).
- **Every planted defect caught** (each `expect_failure()` passed): target 5 (permuted t shrunk by 1e-6: A's
  `#{p <= 1/3}` = 12,696 > 4,608); 4 (i) (`qvalue::qvalue` ignoring pi0) and 4 d-i (BH over the whole family for the
  kept one); R-A1, R-A3, R-B2, R-C1 (ticket's), R-D2, R-F1, R-G0, R-H1, R-J2, R-K1, R-L1, R-M1, R-N1.

**Where a pair departs from `models.md` or the ticket text (mechanics, for the auditor; no expected value changed).**
- R-A3: the ticket's call `computeProtRuvModelFactors(…, ruv_k_rule = "user", ruv_k = 0)` is refused by the code
  (`validateProtRuvModelK`, k >= 1; the existing test `test-prot-ruv-factors-in-model.R:307` expects that refusal).
  The k = 0 branch (`R/func_prot_norm_optimization_helpers.R:1352`) is reached through the estimated count
  (`ruv_k_rule = "num_sv_be"`, `estimateRuvFactorCount` mocked to return 0). The DAs at k = 0 are not re-driven here;
  the existing tests drive them (`test-prot-ruv-factors-in-model.R` "an object without the declaration ...",
  `test-kit127-omics-ruv-factors.R`).
- R-F1's planted defect "`Y - 2 * W alpha`" is read with `W` the library's uncentred W (`RuvIII.W` in `models.md`): with
  the centred W_c the defect stays in col(W_c) and R-F1 itself says it changes no group coefficient.
- R-M1's planted defect: `[X, W Q]` with a rank-1 Q of two columns is not of full rank and edgeR refuses it, so the
  planted design is `[X, W q]`, one column (the same rank-1 map).
- R-B1, R-B2, R-E1 and R-L1 take `alpha` from `computeRuvCovariates()$alpha` (ruv's); the quantities used are
  invariant under a change of basis of U, and RUV-III-C's U spans the same space on complete data (R-J1).

Core-hours: 40.8 job minutes, about 2.72 (job minutes x 4 / 60; set-up included in the job time).

## Resume here

State after the fourth implementer session (2026-11-24): all 142 frozen statements proved; 142 theorems listed in
`axioms.tsv`; audited sound (proofs audit 3; the last three Proof modules were `Proof/RemovedFactors.lean`,
`Proof/RuvIIICPattern.lean`, `Proof/RuvControlRules.lean`). Pairs passed: target 3 (P0a), target 1 (P1a), and in run
45335802998 target 5 (P2), target 4's library parts and group R on merged code (P0b), each with its planted defect
caught. The owed digests job is done: `bindings.tsv` has 6 rows `pending` (target 2's five, KIT-190; T4
`processEnrichments`, KIT-203). Shipped P at (320, 137) = 57188 (`p2-p0b/p-320-137.tsv`): the first counterexample
of 6.1 (g) is live in the shipped arithmetic.

Next, each waiting on what is named:
1. After KIT-203 merges: target 4 (ii) and the `processEnrichments` flag binding (ticket 8.1).
2. After KIT-190 merges: target 2's pair (`pending/test-kit193-pairs-replay.R`, with the builder comparison still to
   write), the L1 cases (the lipid method in R-G0 and L1's control selection), R-N1's sweep-helper cases with
   `getNegCtrlProtAnovaHelper`'s binding row; O1 (a statement for the lipid and metabolite rule) is the
   coordinator's call.
3. After KIT-184 merges: R-I1's pair and its bindings (P3), with ticket 8.2's diff of the merged writer.
4. The Permutation oracle's bridge lemma (above) into `Proof/` and `Audit.lean`, by a session allowed to write proofs.
5. The tripwire `test-kit193-bindings.R` (not written); R CMD check base against head (A7) at P3; the record job L1 at
   the final commit.
Lean calls: `lr.sh`-style wrapper logging `LEAN_NUM_THREADS=4 taskset -c 0-3 timeout 900
/home/doktersmol/.elan/bin/lake …` per call in `cpu-time.tsv`; one at a time.

## Not checked

- R was run only on GitHub runners (P0a, P1a, P2/P0b).
- The proofs are checked on the laptop only (`lake build --wfail`, `Audit.lean`); the record job L1 runs once, at the
  final commit. (Two earlier lines here, on the gate's reads and on the CI workflow, were stale: the statements register
  records the gate's reads, and the workflow and L0 are above.)
- Target 2's pair, target 4 (ii), the L1 cases, R-N1's sweep cases and R-I1's pair have not run (they wait for
  KIT-190, KIT-203 and KIT-184). R CMD check (A7) has not run on this branch. The Permutation oracle's equality to
  the Spec's counts is argued, not yet a Lean theorem.
