### KIT-062 MOFA2 diagnosis: jobs

| job | combination | CPU | numpy BLAS threads | float32 option (default) | iterations | train s | ELBO last | Z1 | Z2 | Z3 | factors md5 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| diag-default-float32-1 | default-float32 | AMD EPYC 9V45 96-Core Processor | openblas:4;openmp:4;openblas:4;openblas:4;openmp:4 | data$use_float32 (True) | 186 | 1.05 | -27000.324394 | 0.9937 | 0.9936 | 0.9892 | 29f363abdabc |
| diag-default-float32-2 | default-float32 | AMD EPYC 7763 64-Core Processor | openmp:4;openblas:4;openmp:4;openblas:4;openblas:4 | data$use_float32 (True) | 186 | 1.83 | -27000.325175 | 0.9937 | 0.9936 | 0.9892 | 1a5340d26d68 |
| diag-default-float32-3 | default-float32 | AMD EPYC 9V45 96-Core Processor | openblas:4;openmp:4;openmp:4;openblas:4;openblas:4 | data$use_float32 (True) | 186 | 0.97 | -27000.324394 | 0.9937 | 0.9936 | 0.9892 | 29f363abdabc |
| diag-default-float32-4 | default-float32 | AMD EPYC 9V74 80-Core Processor | openblas:4;openblas:4;openmp:4;openmp:4;openblas:4 | data$use_float32 (True) | 191 | 1.4 | -26999.648357 | 0.9938 | 0.9938 | 0.9892 | b92a8ad49f27 |
| diag-default-float64-1 | default-float64 | AMD EPYC 9V45 96-Core Processor | openmp:4;openblas:4;openblas:4;openmp:4;openblas:4 | data$use_float32 (True) | 186 | 1.06 | -27000.322787 | 0.9937 | 0.9936 | 0.9892 | d41fd116789a |
| diag-default-float64-2 | default-float64 | AMD EPYC 9V74 80-Core Processor | openmp:4;openblas:4;openmp:4;openblas:4;openblas:4 | data$use_float32 (True) | 186 | 1.36 | -27000.322787 | 0.9937 | 0.9936 | 0.9892 | 7a326abb0c53 |
| diag-default-float64-3 | default-float64 | AMD EPYC 9V74 80-Core Processor | openblas:4;openblas:4;openmp:4;openmp:4;openblas:4 | data$use_float32 (True) | 186 | 1.36 | -27000.322787 | 0.9937 | 0.9936 | 0.9892 | 7a326abb0c53 |
| diag-default-float64-4 | default-float64 | AMD EPYC 9V74 80-Core Processor | openmp:4;openmp:4;openblas:4;openblas:4;openblas:4 | data$use_float32 (True) | 186 | 1.77 | -27000.322787 | 0.9937 | 0.9936 | 0.9892 | 59588a0ca70e |
| diag-single-float32-1 | single-float32 | AMD EPYC 9V74 80-Core Processor | openmp:1;openmp:1;openblas:1;openblas:1;openblas:1 | data$use_float32 (True) | 191 | 1.33 | -26999.648357 | 0.9938 | 0.9938 | 0.9892 | b92a8ad49f27 |
| diag-single-float32-2 | single-float32 | AMD EPYC 9V74 80-Core Processor | openblas:1;openblas:1;openblas:1;openmp:1;openmp:1 | data$use_float32 (True) | 191 | 1.31 | -26999.648357 | 0.9938 | 0.9938 | 0.9892 | b92a8ad49f27 |
| diag-single-float32-3 | single-float32 | AMD EPYC 9V74 80-Core Processor | openblas:1;openmp:1;openmp:1;openblas:1;openblas:1 | data$use_float32 (True) | 191 | 1.24 | -26999.648357 | 0.9938 | 0.9938 | 0.9892 | b92a8ad49f27 |
| diag-single-float32-4 | single-float32 | AMD EPYC 7763 64-Core Processor | openblas:1;openmp:1;openblas:1;openmp:1;openblas:1 | data$use_float32 (True) | 191 | 1.74 | -26999.648643 | 0.9938 | 0.9938 | 0.9892 | 073ec410ff95 |
| diag-single-float64-1 | single-float64 | AMD EPYC 9V74 80-Core Processor | openmp:1;openblas:1;openblas:1;openmp:1;openblas:1 | data$use_float32 (True) | 186 | 1.36 | -27000.322787 | 0.9937 | 0.9936 | 0.9892 | 7a326abb0c53 |
| diag-single-float64-2 | single-float64 | AMD EPYC 9V74 80-Core Processor | openblas:1;openmp:1;openblas:1;openmp:1;openblas:1 | data$use_float32 (True) | 186 | 1.27 | -27000.322787 | 0.9937 | 0.9936 | 0.9892 | 7a326abb0c53 |
| diag-single-float64-3 | single-float64 | INTEL(R) XEON(R) PLATINUM 8573C | openblas:1;openmp:1;openmp:1;openblas:1;openblas:1 | data$use_float32 (True) | 186 | 1.69 | -27000.322787 | 0.9937 | 0.9936 | 0.9892 | 5607e3fb53ac |
| diag-single-float64-4 | single-float64 | Intel(R) Xeon(R) Platinum 8370C CPU @ 2.80GHz | openmp:1;openblas:1;openblas:1;openblas:1;openmp:1 | data$use_float32 (True) | 186 | 1.85 | -27000.322787 | 0.9937 | 0.9936 | 0.9892 | bcd3619d9b20 |

### Within default-float32

| pair | same CPU | bit-identical | iterations | factors | weights | variance explained | ELBO | r2 equal (3 dp) | planted equal (3 dp) | min matched factor abs r |
|---|---|---|---|---|---|---|---|---|---|---|
| diag-default-float32-1 vs diag-default-float32-2 | False | False | 186/186 | 1.91e-06 (1.1e-06 rel) | 1.67e-06 (1.06e-06 rel) | 1.19e-05 (4.2e-07 rel) | 0.000782 (2.89e-08 rel) | True | True | 1.0 |
| diag-default-float32-1 vs diag-default-float32-3 | True | True | 186/186 | identical | identical | identical | identical | True | True | 1.0 |
| diag-default-float32-1 vs diag-default-float32-4 | False | False | 186/191 | 0.0194 (0.0112 rel) | 0.0165 (0.0105 rel) | 0.00758 (0.000267 rel) | 0.676 (2.5e-05 rel) | False | True | 0.999998 |
| diag-default-float32-2 vs diag-default-float32-3 | False | False | 186/186 | 1.91e-06 (1.1e-06 rel) | 1.67e-06 (1.06e-06 rel) | 1.19e-05 (4.2e-07 rel) | 0.000782 (2.89e-08 rel) | True | True | 1.0 |
| diag-default-float32-2 vs diag-default-float32-4 | False | False | 186/191 | 0.0194 (0.0112 rel) | 0.0165 (0.0105 rel) | 0.00758 (0.000267 rel) | 0.677 (2.51e-05 rel) | False | True | 0.999998 |
| diag-default-float32-3 vs diag-default-float32-4 | False | False | 186/191 | 0.0194 (0.0112 rel) | 0.0165 (0.0105 rel) | 0.00758 (0.000267 rel) | 0.676 (2.5e-05 rel) | False | True | 0.999998 |

### Within default-float64

| pair | same CPU | bit-identical | iterations | factors | weights | variance explained | ELBO | r2 equal (3 dp) | planted equal (3 dp) | min matched factor abs r |
|---|---|---|---|---|---|---|---|---|---|---|
| diag-default-float64-1 vs diag-default-float64-2 | False | False | 186/186 | 1.08e-14 (6.18e-15 rel) | 5.55e-15 (3.53e-15 rel) | 3.55e-14 (1.25e-15 rel) | 7.28e-12 (2.69e-16 rel) | True | True | 1.0 |
| diag-default-float64-1 vs diag-default-float64-3 | False | False | 186/186 | 1.08e-14 (6.18e-15 rel) | 5.55e-15 (3.53e-15 rel) | 3.55e-14 (1.25e-15 rel) | 7.28e-12 (2.69e-16 rel) | True | True | 1.0 |
| diag-default-float64-1 vs diag-default-float64-4 | False | False | 186/186 | 2.13e-14 (1.22e-14 rel) | 1.47e-14 (9.32e-15 rel) | 5.68e-14 (2e-15 rel) | identical | True | True | 1.0 |
| diag-default-float64-2 vs diag-default-float64-3 | True | True | 186/186 | identical | identical | identical | identical | True | True | 1.0 |
| diag-default-float64-2 vs diag-default-float64-4 | True | False | 186/186 | 2.84e-14 (1.63e-14 rel) | 1.92e-14 (1.22e-14 rel) | 9.24e-14 (3.25e-15 rel) | 7.28e-12 (2.69e-16 rel) | True | True | 1.0 |
| diag-default-float64-3 vs diag-default-float64-4 | True | False | 186/186 | 2.84e-14 (1.63e-14 rel) | 1.92e-14 (1.22e-14 rel) | 9.24e-14 (3.25e-15 rel) | 7.28e-12 (2.69e-16 rel) | True | True | 1.0 |

### Within single-float32

| pair | same CPU | bit-identical | iterations | factors | weights | variance explained | ELBO | r2 equal (3 dp) | planted equal (3 dp) | min matched factor abs r |
|---|---|---|---|---|---|---|---|---|---|---|
| diag-single-float32-1 vs diag-single-float32-2 | True | True | 191/191 | identical | identical | identical | identical | True | True | 1.0 |
| diag-single-float32-1 vs diag-single-float32-3 | True | True | 191/191 | identical | identical | identical | identical | True | True | 1.0 |
| diag-single-float32-1 vs diag-single-float32-4 | False | False | 191/191 | 2.03e-06 (1.16e-06 rel) | 2.15e-06 (1.37e-06 rel) | 1.19e-05 (4.2e-07 rel) | 0.000287 (1.06e-08 rel) | True | True | 1.0 |
| diag-single-float32-2 vs diag-single-float32-3 | True | True | 191/191 | identical | identical | identical | identical | True | True | 1.0 |
| diag-single-float32-2 vs diag-single-float32-4 | False | False | 191/191 | 2.03e-06 (1.16e-06 rel) | 2.15e-06 (1.37e-06 rel) | 1.19e-05 (4.2e-07 rel) | 0.000287 (1.06e-08 rel) | True | True | 1.0 |
| diag-single-float32-3 vs diag-single-float32-4 | False | False | 191/191 | 2.03e-06 (1.16e-06 rel) | 2.15e-06 (1.37e-06 rel) | 1.19e-05 (4.2e-07 rel) | 0.000287 (1.06e-08 rel) | True | True | 1.0 |

### Within single-float64

| pair | same CPU | bit-identical | iterations | factors | weights | variance explained | ELBO | r2 equal (3 dp) | planted equal (3 dp) | min matched factor abs r |
|---|---|---|---|---|---|---|---|---|---|---|
| diag-single-float64-1 vs diag-single-float64-2 | True | True | 186/186 | identical | identical | identical | identical | True | True | 1.0 |
| diag-single-float64-1 vs diag-single-float64-3 | False | False | 186/186 | 4.33e-15 (2.49e-15 rel) | 4.66e-15 (2.96e-15 rel) | 3.55e-14 (1.25e-15 rel) | identical | True | True | 1.0 |
| diag-single-float64-1 vs diag-single-float64-4 | False | False | 186/186 | 2.2e-14 (1.26e-14 rel) | 1.21e-14 (7.69e-15 rel) | 6.75e-14 (2.38e-15 rel) | 1.09e-11 (4.04e-16 rel) | True | True | 1.0 |
| diag-single-float64-2 vs diag-single-float64-3 | False | False | 186/186 | 4.33e-15 (2.49e-15 rel) | 4.66e-15 (2.96e-15 rel) | 3.55e-14 (1.25e-15 rel) | identical | True | True | 1.0 |
| diag-single-float64-2 vs diag-single-float64-4 | False | False | 186/186 | 2.2e-14 (1.26e-14 rel) | 1.21e-14 (7.69e-15 rel) | 6.75e-14 (2.38e-15 rel) | 1.09e-11 (4.04e-16 rel) | True | True | 1.0 |
| diag-single-float64-3 vs diag-single-float64-4 | False | False | 186/186 | 1.93e-14 (1.11e-14 rel) | 9.88e-15 (6.28e-15 rel) | 6.75e-14 (2.38e-15 rel) | 1.09e-11 (4.04e-16 rel) | True | True | 1.0 |

### Across combinations (first job of each)

| pair | same CPU | bit-identical | iterations | factors | weights | variance explained | ELBO | r2 equal (3 dp) | planted equal (3 dp) | min matched factor abs r |
|---|---|---|---|---|---|---|---|---|---|---|
| diag-default-float32-1 vs diag-default-float64-1 | True | False | 186/186 | 4.32e-06 (2.48e-06 rel) | 2.51e-06 (1.59e-06 rel) | 1.89e-05 (6.67e-07 rel) | 0.00161 (5.95e-08 rel) | False | True | 1.0 |
| diag-default-float32-1 vs diag-single-float32-1 | False | False | 186/191 | 0.0194 (0.0112 rel) | 0.0165 (0.0105 rel) | 0.00758 (0.000267 rel) | 0.676 (2.5e-05 rel) | False | True | 0.999998 |
| diag-default-float32-1 vs diag-single-float64-1 | False | False | 186/186 | 4.32e-06 (2.48e-06 rel) | 2.51e-06 (1.59e-06 rel) | 1.89e-05 (6.67e-07 rel) | 0.00161 (5.95e-08 rel) | False | True | 1.0 |
| diag-default-float64-1 vs diag-single-float32-1 | False | False | 186/191 | 0.0194 (0.0112 rel) | 0.0165 (0.0105 rel) | 0.00759 (0.000268 rel) | 0.674 (2.5e-05 rel) | False | True | 0.999998 |
| diag-default-float64-1 vs diag-single-float64-1 | False | False | 186/186 | 1.08e-14 (6.18e-15 rel) | 5.55e-15 (3.53e-15 rel) | 3.55e-14 (1.25e-15 rel) | 7.28e-12 (2.69e-16 rel) | True | True | 1.0 |
| diag-single-float32-1 vs diag-single-float64-1 | True | False | 191/186 | 0.0194 (0.0111 rel) | 0.0165 (0.0106 rel) | 0.00759 (0.000268 rel) | 0.674 (2.5e-05 rel) | False | True | 0.999998 |

### Against the stored result (package lane run 44459109264)

| job | factors identical | factors max abs diff |
|---|---|---|
| diag-default-float32-1 | False | 0.0194 |
| diag-default-float32-2 | False | 0.0194 |
| diag-default-float32-3 | False | 0.0194 |
| diag-default-float32-4 | True | 0 |
| diag-default-float64-1 | False | 0.0194 |
| diag-default-float64-2 | False | 0.0194 |
| diag-default-float64-3 | False | 0.0194 |
| diag-default-float64-4 | False | 0.0194 |
| diag-single-float32-1 | True | 0 |
| diag-single-float32-2 | True | 0 |
| diag-single-float32-3 | True | 0 |
| diag-single-float32-4 | False | 2.03e-06 |
| diag-single-float64-1 | False | 0.0194 |
| diag-single-float64-2 | False | 0.0194 |
| diag-single-float64-3 | False | 0.0194 |
| diag-single-float64-4 | False | 0.0194 |
