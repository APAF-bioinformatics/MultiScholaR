# CORE-008 outcome: closed without change (premise false)

Ticket: `tickets/openai/CORE-008.md` revision 1.1. Decision 2026-10-08 (user): items 1–2 were closed on measured phase times, item 3 after an in-process profile, and items 4–5 because they are off the measured path and unmeasured. **No code changed; no guarded run spent.**

## Items 1–2 (P14 path)

Phase times from stage `core-010-control` (four members, protocol v7):

| Phase | Seconds | Share of run |
|---|---|---|
| `relational_validation` (item 1) | 1.8–1.9 | about 2% of 85–119 s |
| `relational_aggregate_order_summary_write` (contains item 2's annotation collapse over a few thousand protein groups) | 55–74 | — |
| `staged_result_read_summary_write` | 29–42 | — |

Removing validation entirely would save about 2 s. That is below run-to-run variation: the same control build ran in 88 s and 119 s.

## Item 3 (per-group appends in the peptide QC filters)

`groups-probe.R` (in this directory) profiles `runDIAReferencePeptideGroups` for `peptideMissingValueImputation` on the small fixture scaled to 500 and 2,000 proteins (control library `jEmE2F`, source at HEAD, `Rprof` at 10 ms).

| Proteins (peptide rows) | Step seconds | `appendWorkflowTableBatch` | `appendWorkflowTableWriter` | `next_group` | imputation kernel |
|---|---|---|---|---|---|
| 500 (18,000) | 2.32 | 1.1% | 1.9% | 4.1% | 16.5% |
| 2,000 (72,000) | 4.51 | 1.1% | 3.1% | 5.0% | 33.3% |

The appends that item 3 would buffer are about 1% of the step, or 3% with the writer. The cost that grows with input is the per-group scientific kernel (`dplyr::summarise` per protein group), which is outside this ticket. The filter variant requires identification evidence the fixture lacks, so the probe used imputation, which goes through the same per-group loop and appends.

## Items 4–5 (wide stream, group-index dictionary)

These are on the state-loading and grouped-write paths, not P14. No consumer with a measured cost was identified, and CORE-008's vehicle cannot observe them. They are left unchanged.

## Consequences

- **The per-group kernel cost** in the peptide QC steps is the candidate for future time work. It is recorded for the PERF-107 hand-off (CORE-016), not opened here.
- **CORE-008's lifecycle and phase-timing work** was done by CORE-017.
