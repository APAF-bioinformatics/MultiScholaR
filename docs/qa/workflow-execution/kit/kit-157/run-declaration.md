# KIT-157 bundle runs on run-k157-rc (declared 2026-11-22, before any run)

## Purpose

`fastlimpadpc` rebuilt for x86-64-v3 (fast_limpa_dpc `portable-x86-64-v3` 0ad8c527 = 87a90b89 with
`-C target-cpu=native` replaced by `-C target-cpu=x86-64-v3` on x86_64). Show on GitHub runners that (1) the Linux bundle's
`fastlimpadpc` holds no AVX-512 instruction whatever runner builds it; (2) its DIA_LIMPA values equal the previous pin's under
the fixed rule; (3) the bundle builds, installs, passes the lifecycle and the DIA_LIMPA e2e lane as before (run 44956911834);
(4) a rough time of the previous and the new build on the same CPU.

## Bindings

| | |
|---|---|
| Branch | `run-k157-rc` (run-only), cut from `kit-157-portable-limpa` becf489f (= `origin/performance-v2` 1c7b2a61 + the pin commit) |
| Instrumentation | 41ba192e and 37414840 (run-k156-rc's 141e7b45 and 71669d61 cherry-picked) and this branch's KIT-157 commit: `packaging/verify/k157/` (old-limpa.sh, old-limpa.R, limpa-run.R, limpa-check.R) and the workflow copy. Nothing under R/, tests/, DESCRIPTION, NAMESPACE, src/, inst/ beyond the pin commit |
| Workflow | `KIT-134 distributable`, `workflow_dispatch` on `run-k157-rc` |
| Read token | `DEPS_READ_TOKEN` only in the Linux library install step (as the existing workflow; the previous pin is built in that same step, its output under RUNNER_TEMP, not printed); no cache; no artifact |
| Comparison | each Linux main job: the bundle's R runs `limpa-run.R` once with this pin (bundle library) and once with the previous pin first on the library path; cases: KIT-129 run-path fixture (seed 129, 120 proteins, 5 per group), three S2-size data sets (seeds 1 to 3, 1,000 proteins x 4 precursors, 6 per group), one larger (seed 157, 5,000 proteins, 10 per group); `proteinMissingValueImputationLimpa(engine = "fastlimpadpc")`, one Rayon thread; every slot of the quantified objects compared by `packaging/verify/compare-results.R` (relative 1e-6, absolute floor 1e-9, discrete exact) |

## Runs, ceiling and attempts

1. **Lifecycle run** (instrumentation proof): stage `library`, platform `linux` (one ubuntu-22.04 job: build with the
   previous pin beside it, bundle checks, the KIT-157 step). Pass: `K157OLD ... installed`; both sides print five
   `K157RUN case` lines and a `fastlimpadpc from` line naming different libraries (different md5); the comparison prints its
   `K157CMP` summary; the perturbed copy is read with exactly 1 leaf outside. Ceiling: job 300 min, step 90 min.
2. **Counted run**: stage `full`, platform `both`: ubuntu-22.04 and ubuntu-24.04 main jobs (build, install, lifecycle off,
   hosts, GUI-entry lanes, e2e lanes on 22.04 including `prot_dia_limpa`, the KIT-157 step) and the windows-2025 main job
   (build, installer, install, lifecycle, GUI-entry and e2e lanes without DIA, figure cases). Ceiling: job 300 min.
3. One correction pass (instrumentation faults only). An unfavourable valid result is a finding, kept and not rerun.

## Judgement (fixed now)

- Portability: `NATIVE fastlimpadpc this pin` zmm count 0 on both Linux jobs, with the runner CPU named.
- Values: every double leaf `agree` under the rule; the largest absolute and relative differences reported. A double
  outside the rule stops the work and goes to the user; nothing is re-recorded. A non-double leaf that differs is examined
  and named.
- Bundle: each step's outcome against run 44956911834 (the lifecycle and e2e lanes with the laptop's counts; GUI-entry
  lanes the known performance-v2 families). A new difference is classified (KIT-134 2.0.4).
- Speed: the `K157RUN case` seconds of the two builds in the same job, reported as measured (one run each; rough).

Results below this line after the runs.

## Lifecycle run [44997249408](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44997249408) on cd4189f3 (stage library, linux): comparison proven; one instrumentation fault

ubuntu-22.04, AMD EPYC 9V74 (avx512f yes). `K157OLD fastlimpadpc 87a90b89: installed (19 s)`. Both sides printed five cases
and loaded different libraries (this pin: bundle library, RemoteSha 0ad8c527, md5 aceabd92...; previous: k157-old,
RemoteSha 87a90b89, md5 7a5be652...; RAYON_NUM_THREADS 1). Comparison: 321 leaves, 321 agree, 0 outside, 0 differ;
141 double leaves, 1,264,660 values; largest absolute and relative difference 0 (bit-identical). Perturbed copy: 1 outside,
as required. Seconds (this pin / previous): item5 0.05/0.06, s2_1 0.12/0.12, s2_2 0.11/0.11, s2_3 0.12/0.12, large 0.55/0.53.
Previous pin's library: 247,132 instructions, zmm 1,522, ymm 629, fma 0 (as KIT-156).

**Fault (instrumentation):** the count of this pin's library split the bundle path at its space (`msr moved`) and read
0 instructions. Correction pass (the only one): paths read whole (`find -print0`), and the perturbed summary written to a
file before its first line is printed (the `head` pipe ended R early, harmless). The counted run follows on that commit.

## Results: counted run [44998302733](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44998302733) on a65c8efa (stage full, both)

Read from the job logs (no artifact, no cache). Conclusion "failure" comes only from the three main jobs' GUI-entry lane
steps, as in run 44956911834. Correction pass used once (before this run, for the lifecycle fault); none left.

| | ubuntu-22.04 | ubuntu-24.04 | windows-2025 |
|---|---|---|---|
| runner CPU | AMD EPYC 9V45 (avx512f yes) | Intel Xeon Platinum 8573C (avx512f yes) | - |
| library (fastlimpadpc 0ad8c527), checks, shell, archive/installer, install | passed | passed | passed |
| lifecycle, installed application | passed | passed | passed |
| hosts (M25, M26) | passed | passed | - |
| browser e2e lanes | 88, 190, 16, 42, **86 (proteomics-dia-limpa)**, 98, 153 = laptop | not run (by design) | 88, 190, 16, 153 = laptop (no DIA, D4) |
| GUI-entry lanes | counts identical to run 44956911834 (prot_dia_limpa 34/4, the known KIT-147..153 families) | the same | identical counts to run 44956911834 (lipid norm_log text digest differs in value only, as between runs) |
| figure paths W-D on / off | - | - | passed / passed |

### fastlimpadpc, previous pin against this pin in the same job

| | 22.04 (EPYC 9V45) | 24.04 (Xeon 8573C) | lifecycle 22.04 (EPYC 9V74) |
|---|---|---|---|
| this pin: instructions / zmm / ymm / fma | 246,612 / **0** / 2,730 / 0 | 246,612 / **0** / 2,730 / 0 | (fault, see above) |
| previous pin (native): instructions / zmm / ymm / fma | 250,773 / 1,522 / 629 / 0 | 246,588 / 8 / 2,610 / 0 | 247,132 / 1,522 / 629 / 0 |
| this pin's library md5 | aceabd92... | aceabd92... | aceabd92... |
| values: leaves agree / outside / differ; largest difference | 321 / 0 / 0; 0 | 321 / 0 / 0; 0 | 321 / 0 / 0; 0 |
| perturbed copy (must be 1 outside) | 1 | 1 | 1 |
| seconds, largest case (20,000 precursors x 20 runs), this / previous | 0.39 / 0.38 | 0.53 / 0.52 | 0.55 / 0.53 |

The new build is the same binary on every runner and holds no AVX-512 instruction; the previous pin took the building
CPU's instructions (1,522 zmm on AMD Zen 4/5, 8 on the Intel runner). Values bit-identical (1,264,660 doubles per side).
Times equal within noise (S2-size cases 0.11 to 0.13 s on both builds).

### run-tests 44997354180 (request k157-limpa-becf489f; 55 DIA_LIMPA test files; base 1c7b2a61; check)

Every file `same` against the base; check identical (0 errors, 3 warnings, 3 notes on both). Both sides run the tools of
the head (new pin), so the previous pin's counts come from the full suite run-tests 44967744129 (KIT-155 candidate
e84b8a12, pin 87a90b89): all 55 files have identical passed/failed/errors/skipped/status. The 11 failing files fail
identically there (known state of performance-v2).
