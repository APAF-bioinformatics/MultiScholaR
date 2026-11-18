# KIT-125: RUV designs on public real data, null splits (2026-11-18)

**Result: one contradiction with KIT-124's simulations by the declaration's rule, in DIA_LIMPA; none in the other three omics.**

**The contradiction (a stop by revision 1.2: KIT-129 does not proceed on it; to the user).** DIA_LIMPA, no RUV: KIT-124 scored it "failed" (rule 1 in cell C5, null rate 0.116), and on real data it is within both limits on every split type: random 0.033 with 2 of 100 splits calling, batch-aligned 0.037 with 0 of 100. The declaration foresaw this for the random splits only (the simulations' own global-null cell C4 gave 0.052). It did not foresee it for the batch-aligned splits, which were the real-data counterpart of C5: there too no RUV is within the limits, on a data set with two real culture batches. Whether this stops KIT-129 is the user's call; nothing here decides it.

Everything else is consistent: today's design is outside the limits in all four omics (in DIA_LIMPA in both of its forms, on both split types), and the recommended design is within both limits in all four.

| Omic (data set) | Split type | No RUV: null rate; splits with a call | Today's design | Recommended | Verdict |
|---|---|---|---|---|---|
| Classic proteomics (Zenodo 17607826) | random | 0.042; 4 of 100, within | **0.078; 28**, outside | 0.040; 5, within | consistent |
| DIA_LIMPA (PXD024043) | random | 0.033; 2, within | correction mode **0.083; 25**, outside; covariates mode 0.058; **11**, outside | 0.051; 7, within | **contradiction (no RUV)** |
| DIA_LIMPA | batch-aligned | 0.037; 0, within | correction mode **0.096; 39**, outside; covariates mode 0.070; **20**, outside | 0.053; 6, within | **contradiction (no RUV)** |
| Transcriptomics (Schurch yeast) | random | 0.048; 4, within | **0.147; 36**, outside | 0.056; 9, within | consistent |
| Metabolomics (MTBLS404) | random | 0.053; 12, within | **0.107; 32**, outside | 0.066; 13, within | consistent |
| Metabolomics | batch-aligned | 0.044; 8, within | **0.099; 26**, outside | 0.059; 6, within | consistent |

KIT-124's rules 1 and 7 scored no RUV "passed" in three omics and "failed" in DIA_LIMPA; today's design "failed" and the recommended design "passed" in all four (`declaration.md`).

Rules, data sets, seeds and limits: `declaration.md`, written before any run. Scores: `runs/k125-results.tsv`, made by `k125_aggregate.R` from `runs/<omic>/runs.tsv`.

## What a row means

One public data set per omic. Per split, 12 biological units of one condition are drawn from the pool and put into two false groups of 6; every p under 0.05 and every call is false. 100 distinct splits per split type, one attempt. The limits (revision 1.2): mean null rate at 0.05 at most 0.075; share of splits with any call at most twice no-RUV's on the same splits, or at most 0.10. No confidence bound is put over the splits (they share samples; revision 1.1): the columns "5th to 95th" and "range" describe the 100 per-split null rates, they are not bounds on the mean.

The three designs are KIT-124's simulator arms 1, 3 and 7 (DIA_LIMPA: also arm 4, today's design in its covariates mode) on the simulations' tree (`0b6fcb6d`): no RUV and today's design through the package's kernels, the recommended design through each simulator's own arm-7 code path (its package kernels are not merged).

## Classic proteomics: Zenodo 17607826, 24 of 30 healthy volunteers' plasma (FragPipe MaxLFQ), random splits

| Design | Mean null rate | 5th to 95th | Range | Splits with any call | Within limits | KIT-124, cell C4 (null rate; any call) | KIT-124 rules 1 and 7 | Verdict |
|---|---|---|---|---|---|---|---|---|
| No RUV | 0.042 | 0.010 to 0.146 | 0.004 to 0.234 | 4 of 100 | yes | 0.014; 0.010 | passed | consistent |
| Today's design | **0.078** | 0.020 to 0.184 | 0.005 to 0.399 | **28 of 100** | **no** (both limits) | 0.073; 0.206 | failed | consistent |
| Recommended | 0.040 | 0.014 to 0.085 | 0.010 to 0.132 | 5 of 100 | yes | 0.050; 0.061 | passed | consistent |

About 563 proteins tested per split. Mean k: today 1.3, recommended 1.3 (never 0). No refusal, no crash. Six of the 30 samples are outside the pool by the declared rule (73 to 217 quantified proteins against a median of 552).

## Transcriptomics: Schurch et al. yeast wild type, 42 biological replicates, random splits

| Design | Mean null rate | 5th to 95th | Range | Splits with any call | Within limits | KIT-124, cell C4 | KIT-124 rules 1 and 7 | Verdict |
|---|---|---|---|---|---|---|---|---|
| No RUV | 0.048 | 0.003 to 0.269 | 0.002 to 0.407 | 4 of 100 | yes | 0.010; 0.000 | passed | consistent |
| Today's design | **0.147** | 0.007 to 0.398 | 0.003 to 0.507 | **36 of 100** | **no** (both limits) | 0.056; 0.078 | failed (rule 1 in C5 only) | consistent |
| Recommended | 0.056 | 0.016 to 0.133 | 0.013 to 0.191 | 9 of 100 | yes | 0.056; 0.074 | passed | consistent |

About 6,050 genes tested per split. Mean k: today 2.4, recommended 2.0 (never 0). No refusal, no crash. The case foreseen in the declaration (today's design within the limits on random splits, a contradiction by the letter) did not occur: today's design is well outside them, and worse than in the simulations' global-null cell (0.147 against 0.056). When no RUV calls anything it calls many genes (median 1,168 in its 4 splits); the recommended design's 9 splits with a call have a median of 20 calls.

## Metabolomics: MTBLS404, 184 volunteers' urine, 120 metabolites

| Split type | Design | Mean null rate | 5th to 95th | Range | Splits with any call | Within limits | KIT-124 (C4 null rate; C5 null rate) | KIT-124 rules 1 and 7 | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| random | No RUV | 0.053 | 0.000 to 0.208 | 0.000 to 0.475 | 12 of 100 | yes | 0.012; 0.036 | passed | consistent |
| random | Today's design | **0.107** | 0.008 to 0.301 | 0.000 to 0.575 | **32 of 100** | **no** (both limits) | 0.098; 0.194 | failed | consistent |
| random | Recommended | 0.066 | 0.008 to 0.142 | 0.000 to 0.225 | 13 of 100 | yes | 0.044; 0.054 | passed | consistent |
| batch-aligned (4+2 against 2+4) | No RUV | 0.044 | 0.000 to 0.209 | 0.000 to 0.500 | 8 of 100 | yes | | | |
| batch-aligned | Today's design | **0.099** | 0.017 to 0.310 | 0.000 to 0.483 | **26 of 100** | **no** (both limits) | | | |
| batch-aligned | Recommended | 0.059 | 0.008 to 0.125 | 0.000 to 0.500 | 6 of 100 | yes | | | |

120 metabolites tested per split, none missing. Mean k: today 2.2 and 1.9, recommended 1.8 (0 in 1 and 2 splits). No refusal, no crash. S4 has no valid any-call figure for the global null (its cell C4 is void by its scoring unit), so the any-call column has no simulated counterpart here.

## DIA_LIMPA: PXD024043, 203 untreated single HeLa cells in two culture batches (65 and 138 cells)

| Split type | Design | Mean null rate | 5th to 95th | Range | Splits with any call | Within limits | KIT-124 (C4 null rate; C4 any call; C5 null rate) | KIT-124 rules 1 and 7 | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| random | No RUV | 0.033 | 0.012 to 0.082 | 0.007 to 0.160 | 2 of 100 | yes | 0.052; 0.047; 0.116 | failed (rule 1 in C5 only) | **contradiction** |
| random | Today's design, correction mode | **0.083** | 0.027 to 0.162 | 0.018 to 0.290 | **25 of 100** | **no** (both limits) | 0.190; 0.576; 0.287 | failed | consistent |
| random | Today's design, covariates mode (shipped default) | 0.058 | 0.019 to 0.126 | 0.008 to 0.239 | **11 of 100** | **no** (any-call limit; the null rate is within) | 0.166; 0.564; 0.265 | failed | consistent |
| random | Recommended | 0.051 | 0.025 to 0.089 | 0.018 to 0.140 | 7 of 100 | yes | 0.049; 0.044; 0.051 | passed | consistent |
| batch-aligned (4+2 against 2+4) | No RUV | 0.037 | 0.016 to 0.074 | 0.010 to 0.104 | 0 of 100 | yes | | | **contradiction** |
| batch-aligned | Today's design, correction mode | **0.096** | 0.043 to 0.170 | 0.018 to 0.227 | **39 of 100** | **no** (both limits) | | | consistent |
| batch-aligned | Today's design, covariates mode (shipped default) | 0.070 | 0.029 to 0.125 | 0.017 to 0.182 | **20 of 100** | **no** (any-call limit; the null rate is within) | | | consistent |
| batch-aligned | Recommended | 0.053 | 0.028 to 0.090 | 0.020 to 0.117 | 6 of 100 | yes | | | consistent |

About 1,135 (random) and 1,120 (batch-aligned) proteins tested per split. Mean k: today 1.7 and 1.9, recommended 2.7 and 2.4 (never 0). No refusal, no crash: 800 rows, all `ok`. Run on a GitHub-hosted runner (run 44207066523, attempt 1), not on the laptop; see the run record.

**The contradiction, by the letter of revision 1.2.** No RUV "failed" KIT-124's rule 1 in DIA_LIMPA (in cell C5 only, 0.116) and is within both limits on both split types here, so the rule says "contradiction". For the random splits the declaration foresaw it: the simulations' global-null cell gave 0.052. For the batch-aligned splits it did not: they tie the recorded culture batch to the false groups at a correlation of 1/3, the layout of C5, and no RUV's mean null rate there is 0.037 against the simulation's 0.116, with no split calling anything. The simulation and the real data disagree about no RUV under a batch tied to the groups at this strength. Why was not looked into: the size of the batch difference in this data set was not measured, and C5's simulated factor may be stronger than the real one.

**What the batch-aligned splits show.** Tying the culture batch to the false groups raised today's design in both of its forms and left the other two where they were: correction mode 0.083 to 0.096 and 25 to 39 splits with a call; covariates mode 0.058 to 0.070 and 11 to 20; no RUV 0.033 to 0.037 and 2 to 0; recommended 0.051 to 0.053 and 7 to 6. The direction for today's design is the simulations' (C4 to C5: 0.190 to 0.287, 0.166 to 0.265); the real rates are a quarter to under a half of the simulated ones. The random splits are not free of the batch either: by chance 99 of 100 hold cells of both batches, and in 68 the two false groups hold unequal numbers of the smaller batch's cells.

Other things a reader should not miss here:

- Today's covariates mode, the route's shipped default, is within the null-rate limit on both split types (0.058, 0.070; the simulations gave 0.166 and 0.265). It is outside on the any-call limit alone, and on the random splits by one split (11 of 100 against a limit of 10).
- The recommended design's 7 of 100 on the random splits is above twice no-RUV's 2 of 100 and passes by the 0.10 branch of the limit.

## Read plainly

- **Today's design** is outside the limits on real data in all four omics, as in the simulations: outside both limits in classic proteomics, transcriptomics, metabolomics and DIA_LIMPA's correction mode; outside the any-call limit only in DIA_LIMPA's covariates mode. In transcriptomics it is worse on real data than in the simulations' global-null cell; in DIA_LIMPA it is milder than there (0.083 and 0.058 against 0.190 and 0.166).
- **The recommended design** is within both limits in all four. Its mean null rate is 0.040 (classic proteomics), 0.051 and 0.053 (DIA_LIMPA), 0.056 (transcriptomics) and 0.066 and 0.059 (metabolomics); the simulations gave 0.050, 0.049, 0.056 and 0.044.
- **No RUV** is within the limits in all four, with a mean null rate of 0.03 to 0.05. In three omics that agrees with KIT-124's "passed". In DIA_LIMPA KIT-124 scored it "failed", so by the declaration's rule this is the contradiction stated at the top. In the three other omics the simulations gave 0.01 in the cells with a strong unwanted factor (C4) and 0.05 in the cell without one (C3: 0.050, 0.053, 0.048); the real figures sit at the second.

Things a reader should not miss, none of which crosses a limit:

- Metabolomics, recommended design, random splits: 13 of 100 splits have a call. That is above 0.10 and passes only because it is not above twice no-RUV's 12 of 100. Its mean null rate, 0.066, is the closest to the 0.075 limit of any passing row.
- On real data no RUV itself calls something in 4 to 12 of 100 null splits in classic proteomics, transcriptomics and metabolomics, against 0 to 1% in the simulations' global-null cell (DIA_LIMPA: 2 and 0 of 100 against 4.7%). The cause was not looked into.
- The per-split null rate ranges widely for every design (up to 0.2 to 0.5 in single splits). The limits are read against the mean, as revision 1.2 fixed.
- Metabolomics: the batch-aligned splits of MTBLS404 gave the same picture as the random ones: tying the recorded batch to the false groups at a correlation of 1/3 did not raise any design's rate. So this data set did not put the designs under the condition where the simulations show today's design failing most (cell C5). The deposited table is the study's processed table; whether batch was already corrected in it was not verified at the source.

## Limits

- A null split measures false calls, not power (revision 1.1). Nothing here says the recommended design finds more.
- One data set per omic is one unit; 100 splits of it share samples.
- Kernel level, on the simulations' tree: the unmerged package kernels of the recommended design (`kit-126`, `kit-127`, `kit-129-k-dia`, `kit-129-trx`) are not what ran.
- DIA_LIMPA ran on a GitHub-hosted runner, the other three on the laptop. The runner's lifecycle splits equal the laptop's (run record), so its figures stand beside the others; no memory or timing figure is taken from it.
- The verdict rule reads means over 100 splits that share samples against fixed limits, with no bound; two rows pass or fail by one split (DIA_LIMPA covariates mode, random: 11 against 10; metabolomics recommended, random: 13 against 24 by the doubling branch only).
- Preparation is minimal and declared (no between-sample normalisation beyond what the deposited tables carry; the GUI route's filters and normalisation are not run).
- Metabolomics has 120 features, so a split's null rate moves in steps of 0.008.

## Not checked on real data

Accepted in the survey (revision 1.3) and not run; the four omics above are the declared pass.

| Workflow | Data set accepted in revision 1.3 | State |
|---|---|---|
| DIA with IQ (classic) | PXD024043 | not run in this pass |
| LFQ MaxQuant | PXD002854 (technical replicates) | not run in this pass |
| TMT | Zenodo 14767905 (technical replicates; needs a PSM-to-protein roll-up) | not run in this pass |
| Lipidomics | ST003456 (CC BY-NC-ND) | not run in this pass; the stub asks for lipidomics or metabolomics |

## Differences from the accepted survey

- Classic proteomics: 24 of the 30 healthy volunteers are in the pool (the survey counted 30); six runs have 73 to 217 quantified proteins.
- Transcriptomics: 42 of 48 replicates (the survey named both figures).
- DIA_LIMPA: the report holds 434 runs (the survey read 435 from the stats file); the 203 untreated cells and their two dates are as surveyed. The paper says the two dates are two independent cell culture batches.
- Sizes and licences as surveyed. The Zenodo file's md5 equals the record's.

## Run record

| Omic | Split types | Splits | Wall time | Exit | Refused or crashed rows |
|---|---|---|---|---|---|
| Classic proteomics | random | 100 | 1,724 s | 0 | 0 |
| Transcriptomics | random | 100 | 899 s | 0 | 0 |
| Metabolomics | random, batch | 100 each | 485 s | 0 | 0 |
| DIA_LIMPA | random, batch | 100 each | 4,803 s (hosted runner) | 0 | 0 |

Classic proteomics, transcriptomics, metabolomics: laptop, one R process at a time, `nice -n 19`, one attempt each, 2026-11-18 10:26 to 11:18 (`runs/confirm.log`). Sidecars: `runs/<omic>/meta.txt` (tree commit, libraries, versions, command, seed rule, input sha256). The lifecycle test before them: `lifecycle/` (3 splits per omic and split type; the deliberate failure recorded as a crash and the run completed, in all four omics).

## DIA_LIMPA on the hosted runner: what was checked

GitHub run 44207066523 (workflow `.github/workflows/k125-null-split.yml`, branch `run-k125` at `93ed2fe9`, the tip of `kit-125` then), one job, attempt 1, success, 2026-11-18 00:20 to 02:02 UTC. It was the workflow's first execution, so the points below were checked against the declaration and not taken from the exit status.

| Check | Found |
|---|---|
| Attempt | 1 (`runs/dia_limpa-runner/job.tsv`; the run's own record says the same); script exit 0 after 4,815 s; no rerun |
| Input | the job log's `sha256sum -c` line reads `DIANN1.8_SingleCells_CellCycle.zip: OK`; the sidecar's sha256 `eebe51f3...` is the manifest's |
| Scripts and tree | sidecar: `k125_dia_limpa.R` and `k125_common.R` sha256 equal the committed files; package tree `0b6fcb6d`, no change under `R/`, `tests/`, `src/` |
| Library | R 4.6.1; limma, limpa, ruv, sva (hard checks) and qvalue, statmod, pkgload equal the laptop's; fastlimpadpc, fastduckdb, RUVIIIC at their pinned commits; 11 packages reported as differing in version and 5 test packages absent, none of them hard-checked (`runs/dia_limpa-runner/library.tsv`: arrow, future, gh, httr2, MASS, Matrix, nlme, plotly, RcppParallel, survival, zip). BLAS on the runner is OpenBLAS 0.3.20, one thread |
| Lifecycle splits 9001 to 9003 against the laptop's rows | identical. Under the project's fixed rule (`packaging/verify/compare-results.R` on `performance-v2`: relative 1e-6, absolute 1e-9, discrete exact; the `seconds` column left out as it is not a measure): 18 leaves, 18 agree; the 24 null rates at 0.05 and at 0.01 differ by 0; status, k, features, calls, seeds and both groups' sample names equal (`runs/dia_limpa-runner/lifecycle-compare.tsv`) |
| Splits | 800 rows: splits 1 to 100 of both split types for each of the four arms, all `ok`, none refused or crashed, no missing measure; seeds as declared; 100 distinct partitions per type; every batch-aligned split is 4+2 against 2+4 on the 20201010 batch |
| Artifact | `k125-dia-limpa-1-dia-limpa`, sha256 `e9858e8b8c56d060941bece7e519d85afef188f4066a21ceef9fa971008aa43d`: eleven text files (per-split aggregates, the public run names of each split, sidecars, the script's console lines, the library table). No build output, no credential, nothing of the PRIDE report beyond those |

In the repository: `runs/dia_limpa/` (the run's folder as uploaded), `runs/dia_limpa.out`, and `runs/dia_limpa-runner/` (the job's record, its library table, its lifecycle splits and both comparisons with the laptop's). The runner took 24 s per split against the laptop's 45 s.
