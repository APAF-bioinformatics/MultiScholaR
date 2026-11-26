# Smoke check, send 2 of 3 (declared before it is sent; the first stays on record in this directory and in runs.tsv)

Purpose: the coordinator's rulings on send 1 (defects 1 to 3). An instrument check with no expectation; no prediction of the ticket is judged and
nothing is a result. Binding: these scripts' commit on the `v0.9-staging` head of the day. Decision: proceed to the audit and L1r/L1c, or fix mechanics.

What the re-send must show (otherwise it is a defect list, not a pass):
1. both scripts parse (`parse: both scripts parse`);
2. every E16 prediction met (`E16: 140 checks, 0 differ`: the 11 cases as predicted, including the `tf` case that reads a column of the letters T
   and F as text, the stand-alone `dropped` condition at 20 of 400, and `smokeflag`, where the scorer refuses an input carrying the smoke flag);
3. the six plants and `exitplant` as predicted (E15 rows `itsd`, `itsdpos`, `emptysample`, `ruvskip`, `ruvdrop`, `daerror`, `normplant`, `exitplant`), where
   `ruvdrop` and E3 differ from 35.6 only as this base requires (`dropped` for the restricted assay, E3's import pattern); if `emptysample` or `daerror` crash again, the call
   stack is in `unit-trace-*.txt` and the cause is located;
4. the metabolomics ITSD identity check (`E3_itsd_identity`) and `LAu` (`E10u`) exercised at least once (job B, M12 deciding);
5. `SMOKE_counted_refusal`: a request with K188_SMOKE=1 and no K188_LIFECYCLE exits non-zero before any unit (script 4 of request A; `failed` is its expected status);
6. the smoke flag is in the provenance of every output and in the rows' `smoke` column; checks, timing and provenance exist even for a unit that is killed
   (`checks-*-start.tsv`, `provenance-*.json` at the start, `unit-*-progress.tsv`).

Split: request A (E16 10 min; M12 lp, las and K172c 20 min; L12 l1c 10 min; the refusal check 2 min; job ceiling (42 + 10) x 4 / 60 = 3.47 core-hours) and request B
(M12 deciding alone, 45 min; ceiling (45 + 10) x 4 / 60 = 3.67). Charged to L1r's line. Running total: send 1 about 2.13 core-hours (22 script minutes + 10 set-up);
send 2 expected about 2.0 (A) + 2.7 (B) = 4.7, ceiling 7.14; total expected about 6.8, against the 3.33 the smoke check was declared at (the excess, about 3.5, is
inside Part L's accepted 2,611.55 and is reported to the coordinator, not entered in the ticket).
