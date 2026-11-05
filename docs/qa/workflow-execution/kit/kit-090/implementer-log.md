# Implementer's log (KIT-090)

Every command run and every source read, in order, with one line on why.

1. `ls -la; cat log.md` in the clean-room folder: see which files exist and the log's current state.
2. Read `spec.md`, `test-identity.R`, `load-candidate.R`, `run-tests.R`, `env.sh` (one `cat`): the task, the acceptance test and the helpers.
3. Wrote `func_general_nipals.R` (first draft from the spec, with switchable summation variants) and `explore-impute.R` (black-box output comparison with `mixOmics::impute.nipals`, test generators copied from `test-identity.R`).
4. `source env.sh; Rscript explore-impute.R`: compare 32 summation-order variants' imputed cells with the reference's outputs.
5. `sed` on `func_general_nipals.R`: inline the option helper (load_candidate copies only the three functions); `source env.sh; Rscript explore-impute.R` again.
6. Result of step 5: the variant with every sum as a BLAS product (crossprod / %*%) for numerators, denominators and the normalisation matched all 5 cases bit for bit; every other variant differed in most cells (~1e-14).
7. `sed`: make the BLAS normalisation the default; `source env.sh; Rscript run-tests.R`: first acceptance run.
8. Result of step 7: RESULT tests 8 expectations 77 failed 0 error 0 skipped 0 (PCA fields identical too, with the singular value and var.tot as extended-precision sum()).
9. Read the first 3 lines of `func_general_nipals.R` (tool requirement before rewriting it), then rewrote it as the final version: the switches removed (BLAS everywhere in the iteration), roxygen @noRd and the Wold (1966) / clean-room header added.
10. `source env.sh; Rscript run-tests.R`: acceptance run on the final file.
