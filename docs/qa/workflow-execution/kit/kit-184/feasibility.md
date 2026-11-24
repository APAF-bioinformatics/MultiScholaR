# KIT-184: feasibility gates

## G0 (reading, no run): passed

The locked `ruv` is 0.9.7.1 from CRAN (`packaging/multischolar.lock.json:3627-3637`). Its source tarball (CRAN archive,
`ruv_0.9.7.1.tar.gz`, SHA-256 a0c54e56ba3d8f6ae178ae4d0e417a79295abf5dcb68bbae26c4b874734d98d8) and that of 0.9.7.2
(SHA-256 83556e1a9257b0f2b1e8bbf8eb770ed8224f6f44f8e4308c6296ab4f5c473d23), the version the ticket read, were
unpacked and compared on 2026-11-24 (no R was run):

| File | 0.9.7.1 against 0.9.7.2 | SHA-256 |
|---|---|---|
| `R/RUVIII.R` | byte-identical | 8902176c096dc80ac888b7c81f63d121679adfe7ff9c207704dc95f16ea7a5d8 |
| `R/RUV1.R` | byte-identical | 4119796511e8d14dcabf555f46b458f3216ff6158c6e29d1f03f35169503b4cd |
| `R/residop.R` | byte-identical | d355d05fe8e8f64b74911b4ae34cdebddecd24aacf63613787167d13bc120d3b |

Read in 0.9.7.1: `Y0 = residop(Y, M)` (`RUVIII.R:35`); `fullalpha` from the left singular vectors of `Y0 %*% t(Y0)`
(`:36-37`); `alpha` its first k rows (`:39`); `W = Y[, ctl] %*% t(ac) %*% solve(ac %*% t(ac))` (`:41`);
`newY = Y - W %*% alpha` (`:42`); `RUV1` returns `Y` when `eta` is NULL (`RUV1.R:4-5`); `residop(A, B) = A - B (B'B)^-1
B' A`. The ticket's reading of 2.3 holds for the locked version.

Not checked: that the tarball served by CRAN's archive is the one the lock's library was installed from (the lock
records the version and the source, not a checksum).

## G5 (reading): cannot occur from any stored fixture; not built

The repository holds no integration store (`git ls-files` lists no path under a `.workflow-entities` folder) and no
integration input whose `matrix_role` is `corrected` (the four manifests of `tests/testdata/integration/corpus/` say
`normalised`). `readIntegrationInput()` refuses a contract-v2 input without `matrix_stage` as outdated
(`R/utils_workflow_integration_input.R`, "saved by an earlier development version"). The refusal of ticket 7.4's
seventh row (a store made before KIT-153 whose `views` are corrected matrices) is therefore not implemented.

## G1 (run R1, 45325536775, head 939a9d10, the head side's logs): holds on every route's fixture

The figures are those `expectIntegrationRemovedFactors()` printed (`r1/g1-lines.txt`; one line per data set written by
a test, the count of identical lines first). Each is computed in the test from the two Parquet matrices read back.

| Route (test file) | Data set | Declared k | Written | Rank | Features complete of rows | s_(k+1)/s_1 | Share | Largest relative residual | Status |
|---|---|---|---|---|---|---|---|---|---|
| Proteomics classic, legacy RUV-III-C (`test-workflow-integration-input.R`) | proteomics | 1 | 1 | 1 | 63 of 63 | 3.7e-15 | 1 | 2.97e-09 | exact |
| The same at two factors | proteomics | 2 | 2 | 2 | 63 of 63 | 4.3e-15 | 1 | 2.03e-09 | exact |
| The same after the correlation filter removed a run in which three otherwise complete proteins were missing | proteomics | 2 | 2 | 4 | 63 of 63 | 0.0928 | 0.993248 | 0.805 | approximate |
| Proteomics DIA with limpa, covariates and correction mode | proteomics | 2 | 2 | 2 | 64 of 64 | 4.4e-15 | 1 | 1.39e-13 | exact |
| Lipidomics, legacy RUV-III-C | lipidomics_LCMS_Pos, _Neg | 1 | 1 | 1 | 44 of 44 | 2.6e-15, 2.4e-15 | 1 | 5.4e-13, 3.2e-13 | exact |
| Lipidomics, legacy at two factors | lipidomics_LCMS_Pos, _Neg | 2 | 2 | 2 | 44 of 44 | 6.4e-15, 8.1e-15 | 1 | 4.4e-14, 3.5e-14 | exact |
| Lipidomics, declared design (`groups_all_features`) | lipidomics_LCMS_Pos, _Neg | 1 | 1 | 1 | 43 of 43 | 4.3e-15, 3.9e-15 | 1 | 1.2e-12, 2.6e-13 | exact |
| Metabolomics, legacy RUV-III-C | metabolomics_LCMS_Pos, _GCMS | 1 | 1 | 1 | 44 of 44 | 2.6e-15, 2.4e-15 | 1 | 5.4e-13, 3.2e-13 | exact |
| Transcriptomics, today's grid and the declared design | transcriptomics | 2 | 2 | 2 | 4,622 of 4,622 | 2.3e-15, 2.0e-15 | 1 | 1.1e-13, 7.7e-14 | exact |

Observed:
- On every route's fixture the written factors are `exact`, their number equals the step's k, and s_(k+1)/s_1 is at
  most 8.1e-15, seven orders of magnitude below the bound 1e-8. No stop rule of ticket 11 applies: (a) the declared k
  is the k applied on every route; (b) the removed part is confined to k directions on the complete features; (c) no
  fixture has a corrected value missing where the normalised one is present for an otherwise complete feature (the
  count of complete features equals the record's on every route).
- The largest relative residual off the written factors is below 1e-8 on every route. It is largest on the classic
  proteomics route (2.0e-09 and 3.0e-09), four to five orders above the other routes: inferred, not shown, to be the
  rounding of that route's corrected values (the level of each protein is restored after RUV-III-C). The bound holds
  with a margin of about 3.
- After the correlation filter removed a run (the fixture built for it), the status is `approximate`: rank 4 against a
  declared 2, share 0.993248. This is the case ticket 5.4 describes, on the shipped classic route.
- DIA with limpa: the written factors span the design's `W1..Wk` (asserted, residual below 1e-8).
- Transcriptomics: the written factors span RUV-III's own W recomputed from `fullalpha` (asserted, below 1e-8).
- Lipidomics declared design: the largest sine of the principal angle between the written factors and the DA's stored
  factors (`args$ruv_factors[[assay]]$W`) is 2.8e-12 and 6.9e-12: on this fixture, which has no missing value, the two
  spans coincide. How far they are apart with missing values was not measured (see "Not checked" in `outcome.md`).

## G2: pending run R1b

In R1 no fixture's written input held a feature with a missing run (63 of 63 proteins and 88 of 88 lipids complete),
so G2 (the compiled RUV-III-C gives one W to every complete target while other features have missing values) was not
yet shown on a route. The two fixtures at two factors were given missing values (commit 621a32b5; an instrument fix
under VAL-54: no expected value changed) and run again as R1b. What R1 does show: on the fixture after the correlation
filter the three proteins that were missing in the removed run raise the rank from 2 to 4, so RUV-III-C did correct
them with factors of their own.

## G3 (run L1) and G4

G4: KIT-183's counted outcome is known: O1, `camera` the test of record for every data set (ledger VAL-94; run
45225854502). Cell C is therefore not run, and K-E1 decides in cells A and B. G3: pending (L1).
