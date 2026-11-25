# KIT-240 outcome (branch point; post-rebase runs owed): the integration tabs on by default

Ticket revision 1.3. Branch `kit-240` from `v0.9-staging` 3eeaf100. Every run on `run-tests` (R 4.6.1, reference
BLAS, ubuntu-22.04); run ids, commits and job minutes in `runs.tsv`. Evidence copies under `r1/` and `r2/`.

| Commit | Content |
|---|---|
| 0470164c | tests only: `test-kit240-integration-default.R` (T1, T3, T4, T6, T7, T8, T9); T2 in `test-integration-gui.R` |
| 8cdd1f6e | fix: package default TRUE; `entry.R` and `lifecycle-test.js` (RUV panel only); `missing` in the probe and the Model tab's sentences; marks reworded (values unchanged); `elbo_last`, `elbo_rank`, the Results and record sentences; moved tests; README, NEWS, guide, ledger `VAL-<merge>`; `default-identity.R`, `wording.md` |
| ba6b44d8 | tests only (r1.3): T10, the Enrichment tab's mark |
| fac10b51 | fix (r1.3): `integrationEnrichmentLabel` in `R/mod_integration_enrichment.R`; `mark_provisional` "default not yet confirmed"; the identity script's allowed insertion |

## Observed

- **R1 (A1), 45579469695 at 0470164c.** T1 (option FALSE, disabled tile, the unavailable panel), T3 (the general
  sentence for MOFA2, reticulate and mofapy2), T4 (`missing` NULL; the runner has MOFA2 and reticulate, so the expected
  piece is mofapy2), T6 (a), (b), (d) and T7 (`elbo_last` and `elbo_rank` NULL, no sentences) failed at their base
  values: 36 expectation failures, 0 errors. T6 (c) and T8 passed; T2 (`test-integration-gui.R`, 397) passed; T9
  skipped (no mofapy2 on run-tests). Stop rule 1 did not fire.
- **R2 (A2 at the branch point), 45581871256, 8cdd1f6e against 3eeaf100.** `r2/compare.tsv`: 21 files `same` (among
  them the moved tests in `test-integration-gui.R`, `-record.R`, `-results.R`, `-outputs.R`, each passing on both
  sides with its own assertions), `test-integration-model.R` `count-changed` (263 to 265 passed, 0 failed), the new
  file `only-on-head` (64 passed, 0 failed, T9 skipped). 0 `new-failure`, 0 `message-changed`.
  `test-workflow-checkpoint-harness.R` has 1 error on both sides, identical (stored laptop reference, not this
  ticket's). R CMD check: 0 ERROR, 3 WARNING, 3 NOTE on both; findings only on head: 0 (`r2/check-compare.txt`).
  Stop rule 2 did not fire.
- **Identity script I, in R2** (`r2/default-identity-verdict.txt`, read from the file): all checks TRUE. (i)
  `app_ui(NULL)` head unset = base TRUE byte for byte in both modes; module equal after the printed substitutions
  (settings note 1, Results note 2, mark "(provisional default)" -> "(decided)" 21 per mode, each found in base and
  absent from head); (ii) head FALSE = base unset for both objects and modes, module the unavailable panel; (iii)
  head identical across two processes. This run is at the branch point; the ticket's A2 evidence is the run after the
  rebase (R1-06), still owed.

## Scope addition r1.3 (coordinator's ruling on R1-04)

The Enrichment tab now marks its eight defaults not yet confirmed (correlation, minimum samples, STRING species, KEGG
organism, Reactome organism, smallest and largest pathway, seed) with " (default not yet confirmed)" after the label.
Reason: release users meet these defaults on the tab, and the ticket's "the mark says so" held only in
`integrationSettings()`. Test T10; the identity script allows exactly this insertion (8 per module render in head,
none in base). Run R3 (45587886914; base ba6b44d8 = the second tests-only commit, head fac10b51; evidence under `r3/`):

- T10 failed on base with 11 expectation failures, 0 errors (`r3/failures-base.txt`). The eight labels carried no
  " (default not yet confirmed)" mark: at ba6b44d8 the mark text is still the r1.2 draft, so the test looked for
  "(default awaiting the owner's decision)" after each label and did not find it. The count was 0, not 8. The mark text
  lacked "not yet confirmed" and contained "owner". On head all of T10 passes; the file moved from 65 passed and
  11 failed to 76 passed and 0 failed (`fixed`); T9 skipped on both.
- The four Enrichment-tab files (`test-integration-enrichment.R`, `-enrichment-options.R`, `-enrichment-tests.R`,
  `-removed-factors.R`) are `same`, with 0 failures on either side.
- Identity script against 3eeaf100 (`r3/default-identity-verdict.txt`): all checks TRUE. In each mode the Enrichment
  mark is inserted 8 times on head and never on base. The substitutions are as in R2: settings note 1, Results note 2,
  mark 21.
- 83.5 job-minutes, 5.57 core-hours (declared ceiling 8.00).

## Not checked (so far)

- Runs M and B and the identity run after the rebase (they wait for KIT-234, KIT-239, KIT-228).
- The systemd controller showing the tabs (reading only; release-checklist item).
- A7 (checkpoint full suite), A5 (owner's wording).
