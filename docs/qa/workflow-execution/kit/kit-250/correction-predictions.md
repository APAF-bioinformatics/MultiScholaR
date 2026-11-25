# KIT-250 correction pass: runs declared before they are sent (2026-11-25)

The ticket's one correction pass, following AUDIT-KIT-250-impl-20261125 (I01 (a), I01 (b), I02, I03). Staging 63143aa3 was merged first (8239aae1; ledger: staging's rows kept, KIT-250's row after VAL-259; nothing else conflicted).

## R2' (comparison): tests-only' 43c1a304 against the fix' 77296fcf

- **Purpose:** show each new assertion that can fail on the current head failing there, and passing on the fix; and show no regression.
- **Binding:** base 43c1a304 (the corrected tests on the code of 7b674539 plus the staging merge), head 77296fcf (I02's fix).
- **Files:** the seven files of R2.
- **Settings:** three shards per side, `check: true`.
- **Ceiling:** 600 job-minutes, inside the correction budget of 645.
- **What it decides:** A2 and A3 for the correction pass.

Predictions:

| File | Base (43c1a304) | Head (77296fcf) |
|---|---|---|
| `test-kit250-organism-dialog.R` | 5 failures, 0 errors, listed below | all pass |
| `test-prot-01c-import-module-contracts.R` | 437 / 0 | 437 / 0 (`same`) |
| `test-prot-01d`, `test-module-ci-prot-import`, `test-module-ci-trx-app`, `test-workflow-gui-memory-lifetimes`, `test-workflow-server-shared` | as R2's head | `same` (R2's counts) |
| R CMD check | | findings identical on both sides |

The five expected base failures in `test-kit250-organism-dialog.R`:
- E4: the id is NULL, not "prot_import_organism_cancelled".
- E5: `removedNotes` is `character(0)`.
- E8 step 2: the logged W1 has an empty id.
- E8 step 3: no `removeNotification` call.
- E8 step 5: no `removeNotification` call.

E3's new assertions (static backdrop, keyboard off, no close button) are guards: they pass on both sides.

E8 is new mechanics: the first `testServer` test of this module's observers with a real DIA owner. If it errors on both sides in its mechanics, before any of its judged assertions, that is an instrument defect under VAL-54: fix the mechanics and resend (at most three times). If it errors on the head only, that is a result: stop and report.

## H1 (browser): the corrected head

- **Purpose:** I01 (b), the fix observed in a headless browser on a runner.
- **Binding:** commit under test 77296fcf (the fix'). The records commit adds `h1.R` under this folder only.
- **Settings:** one scripts job, `shinytest2`, test mode, network namespace.
- **Ceiling:** 60 job-minutes (script timeout 45). Expected about 20 job-minutes.
- **Predictions:** written in `h1.R` (`K250_PREDICTIONS`) before the run.
- **What it decides:** whether the fix holds in the application. Any DIFFERS row is a result, reported to the coordinator; it is not rerun.
