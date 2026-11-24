# KIT-194 feasibility gates (ticket section 8)

One request, `kit194-g0-5c727c10` (run 45304509096), on the tests-only commit 5c727c10 (`R/` equal to `v0.9-staging`
447ba678): a `scripts` job running `fj-confirm.R` after `R CMD INSTALL`, and a tests job running
`test-omics-norm-state-guards.R`. R 4.6.1, reference BLAS, ubuntu-22.04. Evidence: `fj-confirm.tsv` (G0), `g0-lines.txt`
(the test file's `KIT-194 F1`, `F2`, `F3`, `effect` and base lines), `g0-tests.tsv`.

## G0: F-j executed on the base (observed)

Every route reproduced on both omics. In each row the run's precondition held (head `<prefix>_normalized`,
`normalization_complete` TRUE), the state the route left at the head has exactly `<prefix>_qc_complete`'s LCMS_Pos matrix
(maximum 1,373,458.9, the pre-normalisation scale; the normalised table's maximum was 20.38), it carries no
normalisation product (products: the intensity filter's, plus the Pearson selection's for J1), the normalisation tab's
status is "complete", and the DA ran on that state and returned 514 result rows.

| Omic | Route | Head before the DA | DA source | DA input max | Normalised max | DA rows | Status |
|---|---|---|---|---|---|---|---|
| lipidomics | J1 (Revert, Apply) | `lipid_correlation_filtered` | same | 1,373,458.9 | 20.38 | 514 | complete |
| lipidomics | J2 (Revert, Skip) | `lipid_norm_complete` | same | 1,373,458.9 | 20.38 | 514 | complete |
| lipidomics | J3 (Apply, Revert, DA) | `lipid_qc_complete` | same | 1,373,458.9 | 20.38 | 514 | complete |
| lipidomics | J4 (Apply, Reset, DA) | `lipid_reset` | same | 1,373,458.9 | 20.38 | 514 | complete |
| metabolomics | J1 | `metab_correlation_filtered` | same | 1,373,458.9 | 20.38 | 514 | complete |
| metabolomics | J2 | `metab_norm_complete` | same | 1,373,458.9 | 20.38 | 514 | complete |
| metabolomics | J3 | `metab_qc_complete` | same | 1,373,458.9 | 20.38 | 514 | complete |
| metabolomics | J4 | `metab_reset` | same | 1,373,458.9 | 20.38 | 514 | complete |

Stop rules: J1 reproduced on both omics (7.4, T6, T7 go on); J3 and J4 reproduced (7.5, T8 go on). Nothing stops.
J5 and J6 are shown on the tests-only side of R1 (T7, T6 (d)); the G0 tests job showed both on the base (`KIT-194 T7
base`: the export names `<prefix>_qc_complete` as its current state; T6 (d): Apply committed from `lipid_normalized`
after a Revert from `lipid_ruv_corrected`).

## F1: the lineage predicate (observed; the predicate written in the test, base has none)

Every case held; no stop. In process (lipid with RUV, metab without): `_itsd_norm` by the products clause; `_log2`,
`_normalized`, `_ruv_corrected`, `_correlation_filtered`, `_norm_complete` by products and dependencies; the intensity
filter, duplicate resolution and Finalize QC made by the wrappers on a normalised head all normalised. Negatives:
`_raw_data_s4`, `_intensity_filtered`, `_qc_complete` before normalisation, `<prefix>_reset` and `trx_reset` all FALSE.
Worker path (the KIT-078 tabs' journey, both omics): the five normalisation states TRUE, the four QC and design states
FALSE. Transcriptomics: `trx_ruvs` TRUE (products), its QC states and `trx_reset` FALSE, the expression filter on
`trx_ruvs` TRUE. The adoption marker is absent on every state the tabs commit (in process and worker). Derivation
(A01): TRUE for Skip at the result, Apply, Apply again, Skip at the correlation state, `_normalized` and `_ruv_corrected`
(in process) and for `_ruv_corrected` and `_correlation_filtered` on the worker path; FALSE for the first run's
`_correlation_filtered` and `_ruv_corrected` after Reset and a second run, and for J1's `_correlation_filtered`.
Labels and the clause that fired per state: `g0-lines.txt`.

## F2: transcriptomics' second legacy-path run (recorded, not a stop)

k 2, then k 1 without Reset, on the base: the second run fails with "Normalization error: duplicate 'row.names' are not
allowed (simpleError)" and `trx_ruvs` keeps `W_1`, `W_2` from the first run. (On the base the Expression Filter on
`trx_ruvs` also committed `trx_expression_filtered` with `W_1`, `W_2` in its design, T2.)

## F3: lipid/metab second run with ITSD by pattern (recorded, not a stop)

No error on either omic; the second run's normalised maximum is 4.33 against the first run's 20.11: ITSD, log2 and the
normalisation applied to the normalised table (normalised twice).

## The test file on the tests-only commit (G0 tests job)

311 passed, 178 failed, 2 errors (`g0-tests.tsv`). F1, T9 and the "both commits" cases of T6 and T8 passed. Three
failures were mechanics of the test file, not results (VAL-54 instrument defects), repaired in c846908d with no expected
value, threshold or stop rule changed; R1's base side (c846908d) is the re-run:
1. T1 selected `<prefix>_normalized` after the `<prefix>_log2` case had truncated the history there (the 2 errors);
   the `_log2` case now comes last.
2. The legacy session files of T1 A08, T3 (f, g) and T8 (b1) held an object read from the owner, which is
   reference-backed and does not load in another store ("invalid hydration reference"); they now hold F172's
   self-contained object.
3. T10's replay comparison counted J1's log2 and normalised states, which a Revert left in the store outside the active
   history; the comparison is now on the active history's states.
