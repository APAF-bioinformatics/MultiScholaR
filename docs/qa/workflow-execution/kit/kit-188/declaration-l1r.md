# L1r: the declaration's base on every cell (DRAFT declaration; not sent; the auditor reads it at the bound commit)

Purpose: on the bound commit, complete every arm and job type on data set 401 of every cell, show E1, E2, E3, E7, E8, E10, E10u, E11, E12 to E14, E3_itsd_identity, the LP check (controls equal E9's positional rule, k the uncapped sweep's), K172c, and measure time (G3). Blind (`K188_LIFECYCLE=1`): estimands written NA; rows' SHA-256 kept. Mechanics counts (samples removed, one-assay removals, statuses) are read for instrument checks only; no arm, condition, threshold, edge or rule changes in response to them without the owner and a ledger row.
Binding: e4a6d4ae plus the scripts' commit (40 digits at send). **L1r's first job runs E16 (179 checks) before any simulation unit; no smoke send exercised the `killed` case.** A miss is an instrument defect of the scorer (VAL-54).
Data sets: 401 for every cell and job type; 402 and 403 for L12 and M12 deciding only (the script has no arm filter). Predictions: the `lifecycle.md` table ("The bound commit") for mechanics fields; E7's keys (`E7_refusal/*`, `E7_invalid_k/*`) are required and read from `predictions-l1c.tsv` (a different invalid-k result is the first finding, not a pass).
Jobs, plan / timeout in minutes (timeout 1.5 times the plan rounded up to 5), core-hours = (minutes + 10 set-up) x 4 / 60:

| Job | Plan / timeout | Expected | Ceiling |
|---|---|---|---|
| L: L12 deciding 401-403; L12, L64 every job type on 401; L24 deciding 401; K172c | 115 / 175 | 8.33 | 12.33 |
| M12: deciding 401-403; reported, lp, las on 401 (E7 unit, `E3_itsd_identity`) | 30 / 45 | 2.67 | 3.67 |
| M168 deciding 401 | 80 / 120 | 6.00 | 8.67 |
| M168 reported, lp, las, BA 401 | 110 / 165 | 8.00 | 11.67 |
| **L1r** | | **25.00** | **36.33** |

Running-ceiling line: used so far 14.72 core-hours (lc2 2.29, smoke series 12.43); L1r's expected 25.00, ceiling 36.33; L1c is charged to its own allowance (declaration-l1c.md). Part L's ceiling and L1r's figures were accepted by the owner (VAL-362, r1.24: Part L 2,640.74); the coordinator re-allocates 2.57 core-hours inside it, so L1r's line reads 51.05 (14.72 used + 36.33 ceiling). Measured: M12 401 units (deciding 1,204 to 1,243 s alone, lp 372, las 596, K172c 229); scaled: all other unit times.
Rule for a unit with no result and no trace: as in `declaration-l1c.md`. The partial-match guard is on in every process; a guard stop is `instrument` and an instrument defect by definition.
