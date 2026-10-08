# Where the P14 grouping memory goes (diagnostic, 2026-10-08)

User-approved diagnostic to choose the order of the remaining memory tickets after CORE-006 and CORE-007 closed without change. It is not a counted run and makes no performance claim.

**Method.** `grouping-memory-probe.R` (in this directory) runs the arm F grouping shape of `prepareDIARelationalAggregate` as a COPY. That shape is: three group keys; ordered lists of raw and normalised values; `list_unique` of modified sequences; three counts. Each variant runs in a fresh R process (DuckDB 1.5.5, eight threads, `preserve_insertion_order`, 12 GiB limit, no spill) on the real `screen$p14` input: `selection-960000`, 960,000 rows, 480,000 groups. The probe reports the growth of peak resident memory during the query. `list_sum` stands in for the package's `dia_r_sum`; both consume the same materialised ordered list. The annotation joins and output ordering of the production query are outside this probe.

| Variant | Query memory (rep 1 / rep 2) | Change | Time (s) |
|---|---|---|---|
| current: ordered value lists, string keys | 458 / 462 MB | — | 1.31 / 1.39 |
| integer keys (CORE-009) | 373 / 370 MB | −19% | 1.17 / 1.20 |
| streaming `sum` instead of ordered lists | 236 / 233 MB | −49% | 0.34 / 0.37 |
| streaming `kahan_sum` (compensated) | 249 / 244 MB | −46% | 0.34 / 0.38 |
| integer keys + streaming `sum` | 152 / 149 MB | −67% | 0.19 / 0.20 |
| `count(DISTINCT)` in place of `list_unique` for forms | 635 / 601 MB | +34% | 1.32 / 1.33 |
| integer keys + streaming `sum` + `count(DISTINCT)` | 151 / 144 MB | −68% | 0.21 / 0.23 |

**Reading.**
- The two ordered value lists account for about half the grouping memory and about three-quarters of its time. This is CORE-012's question: whether an order-free compensated sum meets the predeclared tolerance rule.
- Integer keys account for about a fifth, and combine with a streaming sum to about two-thirds.
- `list_unique` should stay; `count(DISTINCT)` is worse.

**Limits.** One cohort of the small screen population; a full cohort is about 8,000,000 rows. The ranking is the result, not the ratios as a full-size prediction. A streaming sum is not the current exact ordered sum, and any adoption goes through CORE-012's tolerance decision.

**Decision (user, 2026-10-08):** run CORE-012 next, then CORE-009, then CORE-010 and CORE-011; CORE-008, CORE-013 and CORE-014 follow, each premise-checked first.
