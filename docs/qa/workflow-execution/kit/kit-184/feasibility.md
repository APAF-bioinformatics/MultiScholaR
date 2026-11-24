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

## G1, G2 (run R1, head side): pending

The per-route figures printed by `expectIntegrationRemovedFactors()` (status, declared and written k, rank, features,
s_(k+1)/s_1, share, largest relative residual) are copied here from R1's logs.

## G3 (run L1) and G4

G4: KIT-183's counted outcome is known: O1, `camera` the test of record for every data set (ledger VAL-94; run
45225854502). Cell C is therefore not run, and K-E1 decides in cells A and B. G3: pending (L1).
