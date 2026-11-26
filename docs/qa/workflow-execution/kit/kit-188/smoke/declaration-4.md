# Smoke series, send 4 (allowed by the coordinator's ruling on script audit 2, for parse and E16 only; no simulation unit, no journey)

Purpose: parse both scripts and run the extended E16 after the repairs of script audit 2 (GUARD-1 to GUARD-3, E5-1, E16-1, BIND-1, RANGE-1, P13-1). Binding: the scripts' commit
on the `v0.9-staging` head of the day. Nothing is a result and no ticket prediction is judged; decision: proceed to L1c and L1r once their base exists, or hand back. A miss is recorded
and handed back without repair; a further send would need a new ruling.

What it must show:
1. both scripts parse (`parse: both scripts parse`);
2. every E16 case finishes and every prediction is met: `E16: 177 checks, 0 differ` (20 cases, 157 predictions), including the new cases `k2la`, `k1eff`, `nss`, `s5fail`, `excess`, `normstop`, `guardstop`,
   `rangestop`, `e5diff`, and `main`'s `lp` and `las` blocks (the LP-LA and LAs-LA tables and E5 passing);
3. the guard-stop case `guardstop`: the scorer stops on a row of status `instrument` whose class contains `partial-match guard` (the child's and the journey's mapping to `instrument` are by reading only: no journey runs);
4. the stops of `normstop` (norm_ok), `rangestop` (a data set outside the declared range) and `e5diff` (E5 differs) as predicted.

Limit: E16 30 minutes; ceiling (30 + 10) x 4 / 60 = 2.67 core-hours, expected about 1.2. Charged to L1r's line.
