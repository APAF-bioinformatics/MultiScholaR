# CORE-006 outcome: closed without change (premise false)

Ticket: `tickets/openai/CORE-006.md` revision 1.1. Decision 2026-10-08 (user): close with no code change and no counted runs if the aligned-read gain was negligible. **No code changed and no guarded run was spent.**

## Finding 1: the planned change does nothing

The ticket planned to change `ROW_GROUP_SIZE 10240` to `10000` in the staged-result `COPY` so that row groups would line up with the reader's 10,000-row pages. DuckDB 1.5.5 (the pinned engine) treats `ROW_GROUP_SIZE` as a target and closes a row group only at a 2,048-row vector boundary, so a target of 10,000 rounds up to 10,240.

Checked with an ordered 480,000-row `COPY` (Snappy, `preserve_insertion_order = true`) at 1 and 8 threads, for both settings. Every file had 47 row groups of 10,240 rows, with the last at 8,960. The planned one-line change would leave the staged file layout unchanged. The check ran through the `duckdb` R package; production issues the same `COPY` through fastduckdb's connection to the same DuckDB 1.5.5 engine.

This is documented DuckDB behaviour, not a defect in DuckDB or fastduckdb. The ticket's premise came from the review and was not checked against the engine when the ticket was written.

## Finding 2: reading in row-group-aligned pages is not worth doing

The waste the ticket targeted is real: most 10,000-row pages filtered on `.dia_output_order` touch two 10,240-row groups. The alternative fix was to read whole row groups and carry the remainder so that downstream batches stay at exactly 10,000 rows. To bound its gain, `read-benchmark.R` (in this directory) was run in-process, outside the guarded harness and not counted. It uses a synthetic staged file shaped like the P14 stage: the summary columns of the real output, the three aggregates, `.dia_output_order`, 4,000,000 rows (one of a full run's three batches), Snappy, `ROW_GROUP_SIZE 10240`, 391 row groups, 113 MB. Reads use `fastduckdb::read_parquet_fast` with eight threads and the order filter, exactly as `streamDIARelationalResultFile` does. Library: `core_005`.

| Read pattern (one 4M-row file) | rep 1 | rep 2 | rep 3 |
|---|---|---|---|
| Current: 10,000-row pages | 1.87 s | 1.90 s | 1.65 s |
| Aligned 10,240-row reads only (ceiling, no carry-over) | 1.53 s | 1.55 s | 1.54 s |
| Aligned reads with a plain R carry-over to 10,000-row batches (earlier run, reps 1–2) | 6.06 s | 5.74 s | — |

- **The ceiling is about 0.27 s per file**, about 0.8 s per full P14 run. The PERF-108 control members took 87.5 s and 87.0 s, so this is about 1% of the run, the same size as the variation between those two runs. The real gain would be smaller still, because carrying rows across pages costs something. The plain R version cost more than the whole saving.
- **Decision:** no change. Spending the phase timer, a lifecycle re-proof and four measured runs to confirm a gain of about 1% is not justified.

## Consequences

- **CORE-005's pending outcome** is settled by the next measured ticket's control pair on `core_005` (CORE-007), against the PERF-108 `full_p14` receipts, under the same rule (peak ≤ 1.00, elapsed ≤ 1.05).
- **Phase timing and the lifecycle re-proof** move to CORE-008, whose validation-phase metric needs them; the ledger row moves with them.
- **Read-phase cost:** once CORE-008's timing exists, its receipts will show the read phase's real share. A future ticket can be opened only if that share is material.
