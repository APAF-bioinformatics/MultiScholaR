# CORE-011 outcome: closed without change (premise false)

Ticket: `tickets/openai/CORE-011.md` revision 1.1. Decision 2026-10-08 (user): close without code change or counted runs after the premise check below. **No code changed; no guarded run spent.** `overflow.json` was not written, because no overflow member ran.

## What the ticket assumed

Turning on DuckDB spill (an owned `temp_directory` with a bounded size, `memory_limit` unchanged) would let a memory-short run complete more slowly instead of failing.

## What was found

`spill-probe.R` (in this directory) runs the P14 grouping shape on `screen$p14` (`selection-960000`, 960,000 rows) with spill enabled. Each limit runs in a fresh process: DuckDB 1.5.5, eight threads, control library `jEmE2F`.

| Shape | DuckDB `memory_limit` | Result | Process memory growth |
|---|---|---|---|
| ordered value lists (P14) | 4 GiB, 450 MiB | completes, 0.7 s | 463–484 MB |
| ordered value lists (P14) | 400 MiB | completes, 2.2 s | 479–484 MB |
| ordered value lists (P14) | 300, 200, 100, 50 MiB | out-of-memory error | — |
| plain `sum` (contrast) | 100 MiB | completes, 0.16 s | 146–153 MB |

1. **The P14 grouping cannot spill materially.** The ordered value lists that the exact sum requires stay in memory. Spill trims at most about 15%, at three times the time, and fails below about 300 MiB on this input.
2. **`memory_limit` does not bound the process.** At 400 MiB the process still grew by 479–484 MB, so a lower engine limit with spill would not make the job fit a smaller budget.
3. **The normal path is unaffected either way.** The 12 GiB engine limit is above the roughly 7.5 GiB peak, so production spill would never trigger.

## Consequences

- **Smaller budgets** are handled by CORE-010's adaptive complete-sample batches (accepted), not by engine spill.
- **The engine and DBI settings** keep spill disabled. CORE-015's constants item records the 12 GiB engine limit as an admission bound only, not a process bound.
- **No other ticket changes.** The remaining tickets (CORE-008, CORE-013, CORE-014) target time and storage, not peak memory.
