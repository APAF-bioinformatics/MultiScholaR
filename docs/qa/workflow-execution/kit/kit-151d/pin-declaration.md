# KIT-151 d, Q3: does the pin give the same model on an Intel processor? Declaration

Fixed before the counted run. Authority: methods ledger INT-13 and INT-17 (user, 2026-11-21). Branch
`run-mofa-evidence-k151d`; lane `.github/workflows/k151d-evidence.yml`, script `.github/mofa-smoke/k151d.R` (task
`pin`), scorer `.github/mofa-smoke/k151d-score.R` (task `pin`). Script hashes: end of this file. Lifecycle:
run 44785248037 (commit fc5e8042): all three tasks, 8 jobs; every deliberately failing cell was written as failed and its job marked failed, every other cell trained to the end (0 at the round limit), the score job ran on all of them. One scorer wording defect was fixed after it (35c41ca3, pin history line); the evidence script and the lane are unchanged. Q3 lifecycle: job 1 on AMD EPYC 7763, job 99 on Intel Xeon Platinum 8370C (ran fully because Intel); fast cap1000 and hard_full23_default were bit-identical across the two under P and PN, not under U (largest raw difference 4.4e-14, same iterations); every OpenBLAS at 1 thread and Haswell under P and PN; 2 minutes per full job. Evidence only: no package code, default or wording changes.

## Purpose

k151c task 2 (run 44678183859) found the pinned maths library bit-identical across three AMD models (EPYC 7763, 9V74,
9V45) and no Intel processor was assigned, so INT-17 says "an Intel processor not yet tested". The historical
divergence (data seed 151007, N = 20 weak, fast) was between an Intel Xeon 6973P-C and an AMD EPYC 7763. This run
repeats the k151c pin task on enough hosted runners to land on Intel processors, and adds two hard slow cases.

## Bindings and conditions

Software as k151c (MOFA2 1.22.0, mofapy2 0.7.3, numpy 1.26.4, scipy 1.12.0, scikit-learn 1.4.0, R 4.6.1, the
bundle's Python; three OpenBLAS copies loaded: R's system 0.3.20, scipy's, numpy's). Conditions, each a separate R
process with the variables set before R (and so Python) starts, as k151c:

- U, unpinned: nothing set.
- P, pinned (what the installer ships, INT-17): `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
  OPENBLAS_CORETYPE=Haswell`.
- PN, secondary: P plus numpy's dispatched AVX-512 targets that the CPU has switched off (fast case only), to
  separate numpy's own SIMD paths from BLAS if P is not enough.

## Cases

- The k151c case: k151 generator `n20weak@7` (data seed 151007), fast, round limit 1,000, the first run's eight arms
  (none, scale_views, cap1000, cap2000, cap5000, cap10000, weight_views, zscore); under U, P and PN.
- Hard slow case 1, `hard_prev10_none`: k151b generator, scenario prev, data seed 151010, no selection, `scale_views`
  off, slow, round limit 10,000. Planted F2 was 0.086 on an Intel Xeon 6973P-C (k151 follow-up) and 0.874 on an AMD
  EPYC 7763 (k151b), the largest cross-processor difference seen; under U and P.
- Hard slow case 2, `hard_full23_default`: k151b generator, scenario full, data seed 151023, the decided selection
  (INT-12), `scale_views` on, slow, round limit 10,000; F2 near the detection limit (0.70 in k151b); under U and P.

## Launch protocol (declared)

The processor of a hosted runner cannot be chosen. Batches of 20 identical jobs (batch 1: jobs 1 to 20; batch 2:
21 to 40; batch 3: 41 to 60; cap 60). Jobs 1 to 4 always run in full (comparators, normally AMD); every later job
first reads its processor and runs in full only on an Intel one (otherwise it records the processor and stops, about a
minute). After batch 1, if at least one Intel job ran in full, stop; otherwise launch batch 2, then batch 3. If no
Intel processor appears in 60 jobs, say so and stop. Every batch is part of this one counted run; the report gives the
number of jobs it took and every processor seen.

## Measures

Per condition and case: classes of bit-identical models (factors, weights and iteration count) across the full jobs,
with their CPU models; for every pair of jobs, whether the vendors differ, identical or not, iterations, the largest
raw difference of factors and weights, the smallest matched |r|. Instrumentation per job and condition: every cell
trained; under P and PN every loaded OpenBLAS reports 1 thread and kernel Haswell; under PN numpy reports no AVX-512
target on. History (instrumentation, not a gate): pinned fast-case rows against k151c's P rows (any AMD model);
unpinned rows against k151c's U rows and k151b's rows on the same CPU model.

## Readings (declared)

- An instrumentation failure (pin not in force, a cell failed) is reported as such and is not a result.
- P: "bit-identical across Intel and AMD" only if every case has one class across all full jobs and at least one
  Intel and one AMD job ran; "agree to rounding" if iterations agree and the largest raw difference is <= 1e-10;
  otherwise "Intel and AMD give different models". The same for PN.
- U: the Intel-AMD differences as measured (iterations, raw differences, planted |r| by vendor): whether the
  historical divergence recurs.

## What the possible outcomes would mean

- P bit-identical across Intel and AMD: INT-17's "an Intel processor not yet tested" can be replaced by the Intel
  models seen; the pin then gives identical models on every processor model tested (still not a guarantee for
  untested processors).
- P not identical but PN identical: numpy's AVX-512 paths matter on Intel; the shipped pin would need the numpy
  variable too (a decision for the user).
- Neither identical: the pin gives replay on one processor family only; the wording must say so.
- No Intel in 60 jobs: the question stays open; said so.

## Budget

Ceiling 6 core-hours (job wall time x 4 cores, including the jobs that stop early). Expected: a full job about 4 to 6
minutes, an early stop about 1 minute; from the lifecycle (2 minutes for 3 cases), a full counted job about 4 minutes; batch 1 with 4 + about 5 Intel full jobs and 11 early stops about 1.5 to 2 core-hours. One counted run (its batches as above); one correction pass at
most, for an instrumentation failure only.

## Script hashes (sha256, fixed for the counted run)

```
8fe56e9ad458fce63443eafe532e9cd47439ca4de73f020f9ae1885024eb6f7e  .github/mofa-smoke/k151d.R
6261544a4ab84361fcfcb6fc7dbc3b0ee19040a61db5520abc3c47cbbc73ba4e  .github/mofa-smoke/k151d-score.R
04c250d48e264f120042f4e3030388fc308480e62c0d9f2ad13a79c577b0eed4  .github/mofa-smoke/install-sva-library.R
405d84a1b78f7cdf2324ef7d85dc309539e88a1264a1ad2a29ca448c0bea2d17  .github/mofa-smoke/install-ruv-library.R
f4caa249375b225ce1d1cf012b6f3c52d4518ee554e465dfc68159274a1051bd  .github/workflows/k151d-evidence.yml
```
