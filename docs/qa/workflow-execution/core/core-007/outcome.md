# CORE-007 outcome: closed without change (premise false)

Ticket: `tickets/openai/CORE-007.md` revision 1.1. Decision 2026-10-08 (user): close without code change or counted runs after the premise check below. **No code changed; no guarded run spent.** `sites.md` was not written, because no site is being changed.

## What the ticket assumed

About 80 `object.size()` checks walk every string after allocation. This forces lazy (ALTREP) string vectors to materialise, which makes `lazy_strings` pointless and raises R memory. The proposed fix was to replace them with byte estimates, turn on `lazy_strings`, and measure combined peak on `full$p14`.

## What was found

1. **The P14 peak is not in R-side sizing.** PERF-108's full attribution (`perf-108/attribution-and-six-batch-outcome.md`) puts the 7.31 GiB whole-process peak inside the third cohort's DuckDB COPY query: grouping, exact sums, joins, ordering. DuckDB reported a 6.05 GiB buffer peak there. The R-side staged read, summary and write phase that this ticket touches peaked at 4.16 GiB. Even a perfect version cannot lower the combined peak, which is the ticket's primary metric.
2. **`lazy_strings` cannot apply to these tables.** fastduckdb returns lazy vectors only for string columns that are not dictionary-encoded (its `read_parquet_fast` documentation). In the real tables from the CORE-004 lifecycle project, every scientific string column is dictionary-encoded:
   - in the P14 output: `Run`, `Protein.Group`, `Protein.Ids`, `Stripped.Sequence` and the row key (`RLE_DICTIONARY`);
   - in the `screen$p14` input: all except `.multischolar_row_key` and `exclusion_reason` (`PLAIN_DICTIONARY`).
3. **`object.size()` does not materialise anything here.** A 2,000,000-row string column read with `lazy_strings = TRUE` came back as a plain `STRSXP`, not ALTREP. `object.size()` on it took 0.02 s and left the R heap unchanged.

## Consequences

- **CORE-005's pending outcome** is settled by the next ticket that runs a control pair on `core_005`.
- **What's left in this ticket:** the per-iteration sizing of growing containers (`readDIAProteinLineage`, `saveWorkflowGroupIndex`) may cost CPU time, not peak memory. `saveWorkflowGroupIndex` is already covered by CORE-008's group-index item. No separate work is opened.
- **Series order:** this finding, with CORE-006's, shows that the remaining memory gains can come only from the DuckDB grouping phase. The series order is revisited in the plan.

## Correction (2026-10-08, post-hoc audit `CORE-SERIES-POSTHOC-AUDIT-20261008b` P03)

Finding 3 sized a column that came back as a plain `STRSXP` because it was dictionary-encoded, which is finding 2's mechanism. So it shows only that such a column is not lazy. Whether `object.size()` materialises a truly lazy, non-dictionary vector was not tested. Two `screen$p14` input string columns are not dictionary-encoded (`.multischolar_row_key`, `exclusion_reason`) and would be returned lazily. The closure rests on findings 1 and 2.
