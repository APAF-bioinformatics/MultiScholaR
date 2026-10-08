# CORE-012: summation contract decision

Status: **predeclared, not yet measured.** Nothing below this section may be changed once results are visible. A superseding decision needs a new dated entry, never an edit.

## Exact reference

The exact sum of each group is computed by adding the doubles as exact rationals (`gmp::as.bigq`). Relative deviation is measured in 256-bit MPFR (`Rmpfr`). Both come from an isolated library outside the repository (`MultiScholaR-dia-proof-20260911-VhC1n5/core/reference-library`: gmp 0.7.5.1, Rmpfr 1.1.3; system GMP 6.3.0, MPFR 4.2.2). The library is never in `DESCRIPTION` and never on the production path. The self-test (`core-012/reference-selftest.R`) passed all six checks:
- an exact cancellation that R's long-double `sum()` loses;
- the exact rational value of 0.1;
- order independence;
- one-ULP resolution;
- the 2^-50 scale;
- exact zero.

Receipt: [`core-012/reference-receipt.json`](core-012/reference-receipt.json).

## Predeclaration (confirmed by the user 2026-10-08)

Copied verbatim from `tickets/openai/CORE-012.md` revision 1.1, item 2 (SHA-256 of the copied text: `95bdb10782eb8860b9819058aa47800ffec31cec33ec8e3c99615727da5fe921`), before any measurement:

2. **Predeclaration (confirmed by the user 2026-10-08; not to be altered after results are visible).** Written verbatim into `core-012-sum-contract.md` before any measurement:
   - Reference: exact rational sum of each group's values via `Rmpfr` from item 1.
   - Tolerance: per-group relative deviation from the reference at most 2^-50, applied to both the ordered long-double sum and the compensated sum.
   - Downstream decisions recomputed on every fixture and on `full$p14`: intensity cutoff via `ceiling(quantile())`, minimum-peptide filtering, replicate filtering, imputation inputs.
   - Adoption rule: adopt the compensated sum only if every downstream decision is identical on all fixtures and the full corpus, every group is within tolerance, the prototype's combined peak ratio against control is ≤ 0.85, and its elapsed ratio is ≤ 1.05.
   - Any failure on any point retains the ordered sum. A narrow miss on the 0.85 bar after CORE-009 and CORE-010 have reduced the peak is a retain outcome, not grounds for a retry or a changed bar.

**Reading under the 2026-10-08 reorder.** CORE-012 now runs before CORE-009 and CORE-010, so the "narrow miss … after CORE-009 and CORE-010 have reduced the peak" clause has no prior reductions to refer to. The bar is unchanged: the prototype's combined peak ratio against the `core_005` control must be at most 0.85, and elapsed at most 1.05. No threshold, tolerance, reference or decision list is altered by the reorder.

## Order of work (cheapest first)

1. **Fixtures and `screen$p14`, in-process, not counted:**
   - ordered sum against R `sum()`, bit for bit;
   - compensated and ordered sums against the exact reference (per-group maximum relative deviation);
   - downstream decisions recomputed for both sums.

   **Any failure ends the ticket with "retain ordered sum", and no guarded run is spent.**
2. Only if step 1 passes: the prototype branch selected by an explicit harness argument; candidate install; control pair on `core_005` (which also records CORE-005's outcome); candidate pair on `full$p14`; full-corpus downstream decisions.
3. One of exactly two decisions, recorded below with the failing point named on "retain".

## Results

### Step 1: `screen$p14`, in-process (2026-10-08; not counted)

Script `core-012/accuracy-check.R`; raw output `core-012/accuracy-screen.json`. Input: `selection-960000` (960,000 rows), 480,000 groups, every group exactly two values. Ordered sum is R `sum()` on the non-missing values in source order, the C kernel's semantics. Compensated sum is DuckDB 1.5.5 `kahan_sum` at eight threads. The reference is the exact rational sum.

| | Raw | Normalised |
|---|---|---|
| Groups where compensated equals ordered, bit for bit | 464,671 / 480,000 | 464,973 / 480,000 |
| Max relative deviation, ordered sum | 1.11e-16 | 1.11e-16 |
| Max relative deviation, compensated sum | 2.99e-16 | 3.02e-16 |
| Within 2^-50 (8.88e-16) | ordered yes, compensated yes | ordered yes, compensated yes |
| Intensity cutoff `ceiling(quantile())` and retention, 11 percentiles (0.5–99) | identical 11/11 | identical 11/11 |

**Determinism (finding outside the predeclared list).** Repeated at eight threads, `kahan_sum` returned identical values in four runs. At one thread against eight threads, 15,329 raw groups differ: the compensated result depends on how DuckDB combines per-thread partial states. The ordered sum does not depend on threads, because each group is reduced from its source-ordered list. The predeclared rule does not name thread or partition dependence. It is recorded here for the user's decision before any guarded run.


## Decision (user, 2026-10-08)

**Retain the ordered sum.**

The failing point is reproducibility. DuckDB's compensated sum (`kahan_sum`) gives results that depend on the thread count and on how the data is partitioned: 15,329 of 480,000 groups differ between one and eight threads. The ordered long-double sum is independent of both. Scientific output must not depend on execution parameters, because the series verifies outputs by exact hashes, and CORE-010 varies the batch count with available memory. The predeclared tolerance and screen-level decision checks passed, and no threshold, tolerance or rule was altered. The decision rests on a property the predeclaration did not name; the user decided it explicitly after the finding was recorded above.

Not done, because the decision was reached at step 1 (cheapest first): the fixture test file, the prototype branch, the candidate install and the full-corpus pair. No guarded run was spent.

Possible future work, not opened: a deterministic, correctly rounded streaming aggregate (custom C aggregate with an exact accumulator) would remove the value lists without thread dependence. Its per-group state cost is unknown. Consider it only if CORE-009 and CORE-010 leave memory too high.

