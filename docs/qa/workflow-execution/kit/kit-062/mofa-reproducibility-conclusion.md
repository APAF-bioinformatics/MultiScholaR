# KIT-062: MOFA2 reproducibility across machines, conclusion

| Evidence | Result |
|---|---|
| Diagnosis run 44462198417 (16 Linux jobs, five CPU models; `mofa-reproducibility-diagnose.md`) | float32 (MOFA2's default): one of two training paths (186 or 191 iterations) by CPU model and thread count; factors differ by up to 0.019, same solution (matched abs r 0.999998). float64: 186 iterations everywhere, agreement about 1e-14 relative, no time cost |
| Diagnosis run 44466163909 (2 Linux, 2 windows-2025, float64; `mofa-double-windows-diagnose.md`) | 186 iterations, ELBO -27000.322787 on all four; Linux against Windows at most 5.3e-15 relative (factors), 2.8e-15 (weights), 1.3e-15 (variance explained), 4.0e-16 (ELBO); the two Windows jobs bit-identical |
| Package lane run 44465676601 (record, EPYC 7763) | double-precision stored result, 186 iterations; abs r Z1 0.9937, Z2 0.9936, Z3 0.9892 |
| Package lane run 44467071883 (two jobs) | EPYC 9V74: largest relative difference 1.63e-14 (factors), 1.22e-14 (weights), 3.25e-15 (variance explained); EPYC 7763: identical; both pass at the declared 1e-5 |

Decision INT-6 (user, 2026-11-19): the model trains in double precision by default (`precision = "double"`, MOFA2 `data_options$use_float32 = FALSE`); single precision stays selectable, with its consequence stated.
