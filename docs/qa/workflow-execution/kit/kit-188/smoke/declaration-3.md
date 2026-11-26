# Smoke check, send 3 of 3 (the last; declared before it is sent; sends 1 and 2 stay on record in this directory and in runs.tsv)

Purpose and binding as in declaration-2.md: an instrument check with no expectation; nothing is a result and no ticket prediction is judged; scripts' commit merged with the
`v0.9-staging` head of the day. Same split (request A: E16, M12 lp/las/k172c, L12 l1c, the counted-request refusal; request B: M12 deciding alone).
New in send 3: `options(warnPartialMatchDollar = TRUE)` and a global handler that stops on "partial match of" at the top of both scripts (and in run-e16.R), and in `quiet()`.

What send 3 must show (if any item is missed: record and hand back, no further repair; a fourth send needs a new ruling):
1. both scripts parse;
2. ALL E16 cases finish and every prediction is met: `main`, `k1`, `k2phi`, `val248`, `k3`, `k4`, `tf`, `smokeflag`, `itsd`, `gap`, `cellgap`, including the stand-alone `dropped` condition at 20 of 400
   and the corrected `GA not_run` count of 398 (`E16: 140 checks, 0 differ`);
3. `emptysample` and `daerror` as 35.6 predicts for this base (E15 rows met; `emptysample`: status not `instrument`, E12 and E13 pass in both assays; `daerror`: crashed at `da`, RUV columns as `ruvskip`);
4. no partial-match stop in any process (no `partial-match guard` in any console, trace or error);
5. the progress files, `checks-*-start.tsv` and the provenance written for every job;
6. everything send 2 showed, again: `E10`, `E10u`, `E3_itsd_identity`, `E14`, `E1`, `E2`, `E6`, `E8`, `E9`, `E11`, `E15_silence` pass; E3 and `ruvdrop` differ only as this base requires (import pattern; `dropped`); `SMOKE_counted_refusal` exits non-zero before any unit.

Limits: A: E16 25 min, M12 lp/las/k172c 20, l1c 10, refusal 2; job ceiling (57 + 10) x 4 / 60 = 4.47 core-hours; B: 45 min, 3.67. Expected use about 2.1 + 2.1 (the same as send 2, E16 scoring four cases at a time).
