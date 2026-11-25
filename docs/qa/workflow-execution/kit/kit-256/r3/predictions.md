# KIT-256 R3: the repaired G3 (AUDIT-KIT-256-impl-20261126 I01, repair (a)), predictions declared before sending

Purpose: show that G3, now read from a mark set just before its second Assign, fails on the old code and passes on the
fix, so a test guards the direct clear (`R/mod_omic_design_builder.R:667` at 57238dcd).

Binding:
- base: `run-k256-r3-base` fb1b4503 = 37451255 (R/ as v0.9-staging 63143aa3, no fix) with the repaired test file;
- head: `kit-256` 83d6ffd7 = the fix merged with v0.9-staging 42846c0e, with the repaired test file.

One request, the new file only, 1 shard per side, no check. Ceiling 2 x 30 job-minutes = 4.0 core-hours, inside the
ticket's correction budget (705 job-minutes). One run; its result is recorded as it falls (no re-send).

Predictions:
- base: G3 fails by assertion with `actual` NULL (`KIT256-OBS G3: <no message>`), not `c("S1","S2")` as R1 recorded;
  every other row as R1 and R2's base side (E1 x2, E2, E3 x2 fail at their base values; G1, G2, G4 pass); the file
  stays 14 passed, 6 failed, 0 errors.
- head: 20 passed, 0 failed, 0 errors (`KIT256-OBS G3: character(0)`).
- compare verdict: `fixed`.
- Not run (the reason for the change): head with `:667` removed would fail G3.
