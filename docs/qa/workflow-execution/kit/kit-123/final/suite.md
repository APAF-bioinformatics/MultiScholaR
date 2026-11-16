# KIT-123 (RUV k bounded, levels kept): full suite at `652eda20` (2026-11-16)

Baseline: KIT-122 revision 1.1 at `ea236689` ([../../kit-122/final2/suite.md](../../kit-122/final2/suite.md)). Build `browser-smoke.C08ihF` at `652eda20` (R code identical to `5e33c65f`). Between the two: R1 (`9e8db5a7`), R2 (`bee0e705`), the W1 reference tests (`0b6fcb6d`), the re-freeze of two norm test files (`9e9f4f9c`), NEWS, evidence. The KIT-124 simulations were paused (SIGSTOP) for the whole run, so the suite ran on a quiet machine.

`R CMD check --no-tests`: 3 WARNINGs, 3 NOTEs, the same six as the baseline. (Exit 1 is the check-mode aggregation stopping at its header, as in every earlier run.)

| Mode | Files | Passed | Failed | Error blocks | No result |
|---|---|---|---|---|---|
| `load_all`, baseline (with reruns) | 406 | 36,954 | 0 | 2 | 0 |
| `load_all`, this run | 407 | 37,033 | 9 | 2 | 0 |
| check, baseline | 406 | 31,887 | 149 | 237 | 0 |
| check, this run | 407 | 31,966 | 158 | 237 | 0 |

## Every difference (`kit-105/skipcmp.py`, per file)

- **`test-prot-ruv-levels.R`, new file, both modes:** 36 passed, 0 failed (R2's failing-first tests).
- **`test-gui-entry-proteomics.R`, both modes: 9 failures, all in the DIA_LIMPA lane, the first at `ruv_covariates`.** A fresh freeze of the five lanes (`rerun/freeze.log`, `rerun/oracle-diff.txt`): four lanes identical at every step; DIA_LIMPA identical before `ruv_covariates`, and from it on the RUV covariates state, the adjusted matrix and the W design differ, the DA table from `da` on. That is R2 (factors from centred controls), which the 940m and 940r journeys show the same way. The frozen oracle `prot_dia_limpa.json` was missed by the re-freeze `9e9f4f9c`; it is replaced (the other four are not). The file alone afterwards: **193 passed, 0 failed** (`rerun/gui-entry-after-refreeze.log`).
- No other file changes count, status or skips. The two `load_all` error blocks are the diagnostic ones, as before.

## e2e lanes (`e2e.txt`, one attempt each, same build)

All seven pass with the baseline's counts: lipidomics 88, metabolomics 190, transcriptomics 16, DIA 42, DIA_LIMPA 86, enrichment 98, TMT and LFQ 153.

**Decision:** no regression. The one failing file failed at the change KIT-123 declares, and only there; its oracle is re-frozen with that attribution and the file passes.
