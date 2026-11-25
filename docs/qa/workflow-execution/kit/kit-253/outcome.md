# KIT-253 G0: a Reset or Save Design pressed while the normalisation's figures run

Ticket: `tickets/kit/KIT-253.md` r1.2 (Re-check 1 of `AUDIT-KIT-253-plan-20261125`: may start). Commit under test:
v0.9-staging `16112917329a63ab062655f04db5bf1dd7df9729`. This branch (`kit-253-g0`) adds only this folder; no package
code.

## Declared before the run

- **Purpose.** Observe P-a (proteomics DIA, Reset to Pre-Normalization), P-c (proteomics DIA, Save Design) and L-a
  (lipidomics, Reset), each pressed while a named managed action runs, against the ticket's prediction (section 9);
  with two positive controls, P0 and L0 (the same journeys, no press), showing that the instrument sees the target
  completing and a publication.
- **Binding.** The package at 16112917 (the script commit differs from it only in this folder; the script checks
  this and stops otherwise). Driver: KIT-233's `journey.R` read with `git show origin/kit-233-journey:...` (its blob
  hash is recorded in `environment.tsv`). Plain launcher `packaging/electron/entry.R`, process launcher, managed
  worker, test mode, network namespace with loopback only. Fixtures `tests/testdata/e2e/prot_dia` and
  `tests/testdata/e2e/lipid_canonical` (the lipidomics browser lane's settings: ITSD off, normalisation none, RUV
  skipped; proteomics at the rendered defaults, RUV Automatic).
- **Ceiling.** One job, 90 job-minutes (script limit 75 minutes) = 6.0 core-hours; up to three VAL-54 instrument
  reruns, 18.0 core-hours; total 24.0. Expected 45 minutes = 3.0 core-hours.
- **What it decides.** All three press sequences as predicted with both controls proven: stop rule 1, the ticket
  closes with no code and the section 13 candidate row is written. Any press sequence not as predicted (with its press
  inside and its control proven): stop rule 2, stop and report. Target not seen, press not inside the action, or a
  control not proven: stop rule 3, an instrument defect, rerun under VAL-54.
- **Prediction and judging rules**, fixed in `g0.R` (`K253_SEQUENCES`, `K253_NOTICES`, `K253_RULES`) before any value
  is read, and written by the run to `prediction.tsv` and `prediction-rules.txt` before the browser phase.

## Runs

See `runs.tsv`.
