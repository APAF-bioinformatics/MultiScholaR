### KIT-062 MOFA2 diagnosis: jobs

| job | combination | CPU | numpy BLAS threads | float32 option (default) | iterations | train s | ELBO last | Z1 | Z2 | Z3 | factors md5 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| diag-default-float64-1 | default-float64 | AMD EPYC 9V45 96-Core Processor | openblas:4;openblas:4;openmp:4;openblas:4;openmp:4 | data$use_float32 (True) | 186 | 0.91 | -27000.322787 | 0.9937 | 0.9936 | 0.9892 | d41fd116789a |
| diag-default-float64-2 | default-float64 | AMD EPYC 9V74 80-Core Processor | openmp:4;openblas:4;openblas:4;openblas:4;openmp:4 | data$use_float32 (True) | 186 | 1.61 | -27000.322787 | 0.9937 | 0.9936 | 0.9892 | 7a326abb0c53 |
| diag-windows-default-float64-1 | default-float64 | AMD EPYC 9V45 96-Core Processor | openblas:4;openmp:4;openblas:4 | data$use_float32 (True) | 186 | 1.38 | -27000.322787 | 0.9937 | 0.9936 | 0.9892 | ce2141018b00 |
| diag-windows-default-float64-2 | default-float64 | AMD EPYC 9V45 96-Core Processor | openblas:4;openmp:4;openblas:4 | data$use_float32 (True) | 186 | 1.2 | -27000.322787 | 0.9937 | 0.9936 | 0.9892 | ce2141018b00 |

### Within default-float64

| pair | same CPU | bit-identical | iterations | factors | weights | variance explained | ELBO | r2 equal (3 dp) | planted equal (3 dp) | min matched factor abs r |
|---|---|---|---|---|---|---|---|---|---|---|
| diag-default-float64-1 vs diag-default-float64-2 | False | False | 186/186 | 1.08e-14 (6.18e-15 rel) | 5.55e-15 (3.53e-15 rel) | 3.55e-14 (1.25e-15 rel) | 7.28e-12 (2.69e-16 rel) | True | True | 1.0 |
| diag-default-float64-1 vs diag-windows-default-float64-1 | True | False | 186/186 | 2e-15 (1.15e-15 rel) | 3.33e-15 (2.12e-15 rel) | 2.13e-14 (7.51e-16 rel) | 3.64e-12 (1.35e-16 rel) | True | True | 1.0 |
| diag-default-float64-1 vs diag-windows-default-float64-2 | True | False | 186/186 | 2e-15 (1.15e-15 rel) | 3.33e-15 (2.12e-15 rel) | 2.13e-14 (7.51e-16 rel) | 3.64e-12 (1.35e-16 rel) | True | True | 1.0 |
| diag-default-float64-2 vs diag-windows-default-float64-1 | False | False | 186/186 | 9.21e-15 (5.29e-15 rel) | 4.33e-15 (2.75e-15 rel) | 3.55e-14 (1.25e-15 rel) | 1.09e-11 (4.04e-16 rel) | True | True | 1.0 |
| diag-default-float64-2 vs diag-windows-default-float64-2 | False | False | 186/186 | 9.21e-15 (5.29e-15 rel) | 4.33e-15 (2.75e-15 rel) | 3.55e-14 (1.25e-15 rel) | 1.09e-11 (4.04e-16 rel) | True | True | 1.0 |
| diag-windows-default-float64-1 vs diag-windows-default-float64-2 | True | True | 186/186 | identical | identical | identical | identical | True | True | 1.0 |

### Across combinations (first job of each)

| pair | same CPU | bit-identical | iterations | factors | weights | variance explained | ELBO | r2 equal (3 dp) | planted equal (3 dp) | min matched factor abs r |
|---|---|---|---|---|---|---|---|---|---|---|

### Against the stored result (package lane run 44459109264)

| job | factors identical | factors max abs diff |
|---|---|---|
| diag-default-float64-1 | False | 0.0194 |
| diag-default-float64-2 | False | 0.0194 |
| diag-windows-default-float64-1 | False | 0.0194 |
| diag-windows-default-float64-2 | False | 0.0194 |
