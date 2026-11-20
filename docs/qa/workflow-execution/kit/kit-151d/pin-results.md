# KIT-151 d, Q3 results: the pin on Intel processors, scored once

Declaration `pin-declaration.md` (292b8a2d, pushed alone). Lifecycle 44785248037 (fc5e8042). Counted run 44786874024
(9b286f5d), batch 1 only: 20 jobs, 8 on Intel, so by the declared protocol batches 2 and 3 were not launched. 12 jobs
ran in full (jobs 1 to 4 and the 8 Intel jobs), 8 AMD jobs stopped after reading the processor. Job wall time 2,426 s:
**2.7 core-hours** (ceiling 6). Tables: `pin-*.tsv`; scorer output `pin-score.md` (unchanged).

## Processors

| vendor | model | jobs (full) |
|---|---|---|
| Intel | Xeon Platinum 8573C | 7 (7) |
| Intel | Xeon 6973P-C | 1 (1), the model of the historical divergence |
| AMD | EPYC 9V74 | 6 (3) |
| AMD | EPYC 9V45 | 5 (1) |
| AMD | EPYC 7763 | 1 (0) |

## Results (10 cases: k151c's eight fast arms on data seed 151007, two hard slow cases)

| condition | cases with one bit-class across all 12 full jobs | Intel-AMD pairs identical | Intel-AMD pairs with different iterations | lowest Intel-AMD matched abs r |
|---|---|---|---|---|
| U unpinned | 0 of 10 | 0 of 320 | 20 | 0.0006 |
| P pinned (shipped) | 10 of 10 | 320 of 320 | 0 | 1 |
| PN pinned, numpy AVX-512 off | 8 of 8 | 256 of 256 | 0 | 1 |

- **Pinned: bit-identical** (factors, weights and iteration counts) across both Intel models and both AMD models
  that ran, in every case, including the historical case and both hard slow cases.
- **Unpinned:** four bit-classes per case, one per processor model. The Xeon Platinum 8573C and the two AMD models
  agree to rounding (same iterations, raw differences <= 3.8e-13). The Xeon 6973P-C gives a **different model**: other
  iteration counts in 5 of 10 cases and matched factors down to abs r 0.0006 (fast cap5000) and 0.0006 (hard prev@10).
  It reproduces the history: hard prev@10 planted F2 0.086 there against 0.874 on every other processor (k151 follow-up
  on this model: 0.086; k151b on EPYC 7763: 0.874); fast weight_views 166 iterations there against 106 (the first
  k151 run's Intel row: 166). The decided default case (hard full@23, selection 1,000, slow) kept the same iterations and
  planted abs r everywhere (factors in another order on the 6973P-C).
- History check: pinned fast rows 96 of 96 equal to k151c's pinned rows (k151c ran on AMD only).

## Instrumentation

All cells trained; under P and PN every loaded OpenBLAS reported 1 thread and kernel Haswell in every job. **PN
instrumentation failure in 3 jobs** (two 8573C, the 6973P-C): numpy still reported its `AVX512_SPR` target on, so PN
there is not "AVX-512 off"; the PN rows of those jobs are not a test of that condition. This does not touch P, which is
the shipped condition and does not set the numpy variable. No correction pass was used.

## What is shown and what is not

Shown: with OpenBLAS on one thread and the Haswell kernel, MOFA2 training gave bit-identical models on the four processor
models of the counted run (Intel Xeon Platinum 8573C, Intel Xeon 6973P-C, AMD EPYC 9V74, 9V45), and in the lifecycle
on two cases between an Intel Xeon Platinum 8370C and an AMD EPYC 7763 (k151c: EPYC 7763, 9V74, 9V45), and without the pin one Intel model (6973P-C) gives different models on cases near the
detection limit. Not shown: other processors (no Windows, no Apple, no other Intel generations), other library
versions, or that a pinned model is closer to the truth than an unpinned one (the pin fixes the code path, not the
stability of weak factors).

Methods sentence the result supports: "With the numerical library restricted to one thread and a fixed code path
(as shipped), training gave bit-identical models on the six server processor models tested (three Intel, three AMD; one Intel model on two cases only);
without this restriction, one processor model gave a different solution for factors near the detection limit."
