# CORE-009 outcome: closed without change (premise false)

Ticket: `tickets/openai/CORE-009.md` revision 1.1. Decision 2026-10-08 (user): close without code change or counted runs after the premise check below. **No code changed; no guarded run spent.** `peptidoform-semantics.json` was not written, because no encoding change is being made.

## What the ticket assumed

The grouping diagnostic (`../diagnostics/grouping-memory.md`) measured integer keys at −19% of the P14 grouping memory. The ticket would encode sample, protein, peptide and modified keys to integer codes, aggregate on the codes and decode labels on the output.

## What was found

The diagnostic's integer-key variant built the dictionaries and the encoded table as a temporary table **before** its "before" reading of peak resident memory (`VmHWM`). The encoding cost was outside the measured span, and because `VmHWM` is a high-water mark, that earlier step could also hide part of the query's own growth. Production would encode inside the same cohort query.

`end-to-end-probe.R` (this directory) measures both routes with dictionary build, encoding join, aggregation, label decode and ordered COPY all inside the measured span. Same input as the diagnostic (`screen$p14`, `selection-960000`, 960,000 rows, 480,000 groups, payload SHA-256 `e1c1cc32…b9f60d`), control library `jEmE2F`, DuckDB 1.5.5, 8 threads, 12 GiB limit, no spill, fresh process per run. `list_sum` stands in for `dia_r_sum`, as in the diagnostic.

| Route | Query memory growth (3 runs) | Time (s, 3 runs) |
|---|---|---|
| current: string keys | 491 / 497 / 490 MB | 0.83 / 0.92 / 1.26 |
| integer keys, encoding and decode included | 467 / 466 / 466 MB | 0.92 / 0.97 / 0.95 |

Outputs were identical (480,000 rows, `identical()` after renaming the decoded columns).

- **Memory:** about −5% of the grouping step, consistent across runs. On the whole-process P14 peak this is within normal run-to-run variation, so a counted pair could not attribute it.
- **Time:** no clear difference; the string route's spread covers the integer route.

The ticket's expected gain does not survive counting the encoding cost, so the change would not justify its added query complexity or four guarded runs.

## Consequences

- **CORE-005's pending outcome** moves to CORE-010's control pair on `core_005`.
- **CORE-010** recalibrates its estimate against the current string-key layout (`core_005`), not a coded layout.
- **Diagnostic:** a dated correction is added to `../diagnostics/grouping-memory.md`. Its other rows (ordered lists, `count(DISTINCT)`) built the same string table for both sides and are not affected by this flaw.
