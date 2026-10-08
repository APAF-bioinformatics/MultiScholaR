# CORE-012 outcome: retain the ordered sum

Decision (user, 2026-10-08), recorded in [`../core-012-sum-contract.md`](../core-012-sum-contract.md).

- **Reference:** isolated gmp/Rmpfr library; self-test 6/6 ([receipt](reference-receipt.json)).
- **Predeclaration:** frozen and hashed before measurement.
- **Step 1 (`screen$p14`, in-process):**
  - both sums within 2^-50 (ordered ≤ 1.11e-16, compensated ≤ 3.02e-16);
  - intensity cutoff and retention identical at 11 percentiles;
  - compensated equals ordered in 96.8% of groups ([data](accuracy-screen.json)).
- **Deciding finding:** `kahan_sum` depends on the thread count (15,329 groups differ between one and eight threads); the ordered sum does not.
- **Decision:** retain. No production change, no prototype, no install, no guarded run.

Hand-offs:
- **CORE-015:** remove the non-exact DuckDB `sum()` diagnostic arm, which was kept only until this ruling. No prototype branch exists to remove.
- **Rollup shape:** CORE-009 and CORE-010 proceed on the ordered sum. CORE-010's rollup shape becomes final with CORE-010.
- **CORE-005's pending outcome:** moves to CORE-009's control pair on `core_005`.
- **Reference library:** `MultiScholaR-dia-proof-20260911-VhC1n5/core/reference-library` stays outside the repository and may be deleted.
