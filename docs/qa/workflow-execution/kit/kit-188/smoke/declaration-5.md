# Prepared declaration: ONE further smoke send for the guard repair (NOT sent; it needs the coordinator's word)

Purpose: show that the partial-match guard's stop reaches a handler at journey level and at unit level (it was raised from a global calling handler, which no `tryCatch` can catch), and that
the parent records an instrument unit row for a child that left no result. An instrument check with no expectation about any estimand; nothing is a result. Binding: the scripts' commit
on the `v0.9-staging` head of the day. Limited to: parse; E16; the job `guard` on the small metabolomics cell M12 data set 401: three units (`guardjourney`, `guardunit`, `guardcontrol`), each one ordinary journey (arm NO, N0, mask none, about a minute).

What it must show:
1. both scripts parse;
2. `E16: 178 checks, 0 differ` (21 cases, 157 predictions), including the new case `guardunit` (a unit row of status `instrument` whose class carries the guard's text: the scorer stops);
3. `guardjourney`: a journey of status `instrument`, class `partial_match_guard: ...`, the unit completing (exit 0, rows written), check `GUARD_planted` TRUE;
4. `guardunit`: the child ends outside any journey; its unit row is status `instrument` with the text `partial-match guard`, the unit's trace in `unit-...-trace.txt` or the unit result, check `GUARD_planted` TRUE, and `unit_completed` TRUE for it by design;
5. `guardcontrol`: an ordinary journey, status `ok`, `GUARD_planted` TRUE, no guard text anywhere;
6. no other process died without a row, a check or a trace.
Limit: one request, scripts E16 25 min and the guard job 10 min; ceiling (35 + 10) x 4 / 60 = 3.0 core-hours, expected about 1.2.
