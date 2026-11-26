# KIT-262 M3: lifecycle L3r declaration (corrected re-send of L3; committed alone)

**What this is.** L3 (run 45899742324, `m3-lifecycle-declaration.md` at 48d31d4c) did not meet its declaration (`l3/`, `runs.tsv`). The coordinator ruled both failures instrument defects under VAL-54 (a blind rehearsal; no result produced or judged; no expected value, margin or rule changed): **re-send 1 of at most 3**, charged to the same 4.00 core-hour L3 line (about 1.3 used by the first send). Everything in `m3-lifecycle-declaration.md` stands unchanged: purpose, binding (a3c69901 plus `docs/qa/workflow-execution/kit/kit-262/` only), the five scripts, their environments (`OMP_NUM_THREADS=1` in every one), seeds, cells, replicates, R values, the stress floor order `after`, every prediction in its section 3, the pass rule, the budget and the timeouts (5, 20, 55, 10, 30 script-minutes). Only the output directory names change (`k262-l3r/s1` to `s5`).

## The two defects (found in L3's logs and the code)

1. **`score-m3-k262.R`, `tangoBound`** evaluated `f` on the length-2 vector of limits (`all(is.finite(f(lims)))`), so Tango's `if (score(hi) > 0)` received a vector and stopped ("the condition has length > 1"). This stopped S2 before any case ran and failed `FULL_SCORER_ON_ROWS` in S3 and S5.
2. **`m3-k262.R`, `m3Endpoints`** read `ruv_applied`, `ruv_not_applied_reason`, `best_k`, `k_sweep` and `k_cap` from `run$ruv`, but `runWorkflowOmicsReferenceRuv` returns the step's value under `run$ruv$value` and lifts only `optimisation` to the top level (`utils_workflow_omics_declarations.R:2624-2629`). The reads gave NULL, silently NA or FALSE; hence E1 4 of 36, E9 0 of 16 and F2 failing exactly where RUV applied, and C1 always "differing" (so no unit was ever "reused"). This was the declared OBSERVE item "the step-value field names".

## The three edits (and nothing else)

1. `tangoBound` evaluates `f` at each limit separately. New self-test case `tangoBound_scalar_on_known_table`: for the table b = 10, c = 5, n = 100 (difference 0.05) each bound is one finite number that solves Tango's equation (Z(m) = -1.6449 at the upper bound, +1.6449 at the lower, the definition of the inverted test), the lower below and the upper above 0.05. The expected values are fixed by that definition and, beside it, by the existing Python hand values (`bound_upper.b10_c5_n100` 0.118596, `bound_lower.b10_c5_n100` -0.014756).
2. `m3Endpoints` reads those five fields from `run$ruv$value`; `optimisation` stays at the top level. A missing field, a field without the assay `LCMS_Pos`, a missing `value` list or a missing `optimisation` stops with a named message ("M3 instrument: the RUV step's value has no field '...'"), which the unit turns into an instrument row (class `script`), never a silent NULL. `m3Pick` accepts a list or a named vector.
3. S4 prints, for each planted unit, its tag, r, the expected status and class, the observed status and class and `ok` or `MISMATCH` (aggregate, no data value). **No expected class is changed.** If the re-run shows a planted unit's expected class was wrong, that is a changed expectation and L3r stops and goes to the coordinator.

No check is removed, loosened or wrapped; no seed, cell, R value, floor order or prediction changes. E0 (41 of 41), bracket balance and the empty package diff were re-run locally (`m3-local-checks.txt`); no R was run.

## E9, E11 and F2: what stays unobserved

E9 (16 units) and F2 (eight cells and the stress cells) are **not yet observed**: L3's failures there were the accessor defect, not evidence about the product. **E11 ran 0 units in L3 as a consequence of defect 2, not by design and not as a third defect:** E11 runs only where C1 re-uses P's row (`c_source` "reused"), and with `ruv_applied` read as FALSE every unit with a C1 winner was classed "differs" while the others had none. With the fields read correctly, and `K262_M3_FORCE_ARMS=1`, E11 runs on every unit where C1's choice equals P's. It stays "not yet observed" until L3r. The seconds measured in L3 are provisional for the same reason (C1 ran its Manual journey in every unit with a winner).

## Script SHA-256 (at this declaring commit)

- `m3-k262.R` ef6428196c907c932e2dd8b70ad1fb0834a62b2272e92498aacf462677a071eb
- `rules-k262.R` 9547d350211350615a4d6f6e8f992be2307bc753e24790dfdb6cd5a1578d9d0a
- `score-m3-k262.R` c7edbc35461d949132d9618e677fff57e549de915813e1c62cd72004a3646fc3
- `reestimate-m3-k262.R` b6d18e54c08d85480d984b40357c6256acafc4e2373b5af51fbabd2e1972e94d
- `e5hand-m3-expected.tsv` 3e78da07ee78eab1cace85636966ad4629d65aff304a945d914368d1479a13bf

## Running-ceiling line

Used 87.24 of 250.32 plus about 1.3 for the first L3 send = about 88.6; L3r at the L3 line's remaining 2.7 (expected 2.57, 3.91 at spread 1.93 against the 4.00 line: any excess is charged to the instrument reserve, the second of at most three instrument re-sends). No counted job and no counted stage 1.
