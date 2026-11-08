# KIT-110 entry-point suite

## Simulation Setup (run 2): trx_screen alone on c1079dc7

`kit-105/lifecycle/lifecycle.R` (the repaired version, `be0a2601`), R 4.6.1, `load_all`, strict skips, `R_LIBS` kit_105. Log: [lifecycle/trx_screen.txt](lifecycle/trx_screen.txt).

- **Result:** identical to the oracle (flag TRUE, from the ListReporter), 32 s.
- **Traced calls, each as the ticket expected:** `mod_omic_norm_server` 1, `runTrxNormalisation` 2 (the refusal and the run), `runTrxNormInProcess` 2, `runTrxNormPublicationsInProcess` 1, `runWorkflowOmicsPublicationInProcess` 8. Of the 8, the norm step accounts for 5 and the DA and RUVIII steps for 3, as W2 measured before the frame move.

The first attempt ran cleanly, so there was no instrumentation failure.

## Run 3: the trx, lipid and metab entry files in parallel on c1079dc7

| File | Lanes | Pass | Fail | Error | Skip | Elapsed |
|---|---|---|---|---|---|---|
| test-gui-entry-transcriptomics.R | trx_screen | 29 | 0 | 0 | 0 | 44 s |
| test-gui-entry-lipidomics.R | lipid_canonical | 35 | 0 | 0 | 0 | 74 s |
| test-gui-entry-metabolomics.R | metab_lc, metab_gc, metab_combined | 105 | 0 | 0 | 0 | 162 s |

**Decision: the 5 lanes are identical** to their oracles. No oracle changed in KIT-110: `git diff d5244ffe --stat -- tests/testdata/gui-entry` is "". The proteomics lanes do not reach this code (design 6.2), so they were not run.
