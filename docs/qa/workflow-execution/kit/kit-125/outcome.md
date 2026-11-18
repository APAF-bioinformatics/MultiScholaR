# KIT-125: RUV designs on public real data, null splits (2026-11-18)

**Result: no contradiction with KIT-124's simulations in the three omics that ran** (classic proteomics, transcriptomics, metabolomics). In each, no RUV and the recommended design are within the limits of revision 1.2, and today's design is outside them. **DIA_LIMPA has not run:** it was sized above the laptop limit and is prepared for a hosted runner; its row is open. Nothing here is the user's decision.

Rules, data sets, seeds and limits: `declaration.md`, written before any run. Scores: `runs/k125-results.tsv`, made by `k125_aggregate.R` from `runs/<omic>/runs.tsv`.

## What a row means

One public data set per omic. Per split, 12 biological units of one condition are drawn from the pool and put into two false groups of 6; every p under 0.05 and every call is false. 100 distinct splits per split type, one attempt. The limits (revision 1.2): mean null rate at 0.05 at most 0.075; share of splits with any call at most twice no-RUV's on the same splits, or at most 0.10. No confidence bound is put over the splits (they share samples; revision 1.1): the columns "5th to 95th" and "range" describe the 100 per-split null rates, they are not bounds on the mean.

The three designs are KIT-124's simulator arms 1, 3 and 7 on the simulations' tree (`0b6fcb6d`): no RUV and today's design through the package's kernels, the recommended design through each simulator's own arm-7 code path (its package kernels are not merged).

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

## DIA_LIMPA: PXD024043, 203 untreated single HeLa cells in two culture batches: not run

Sized at 2.5 core-hours for the two split types (45 s per split, 40 of them the sweep of today's design), above the declared laptop limit of one core-hour. Prepared for a GitHub-hosted runner; see "To launch" below. The lifecycle test's six splits ran and are in `lifecycle/dia_limpa-reference/`; six splits are not a result and no figure from them is quoted.

## Read plainly

- **Today's design** is outside both limits on real data in all three omics that ran, as in the simulations. In transcriptomics it is worse on real data than in the simulations' global-null cell.
- **The recommended design** is within both limits in all three. Its mean null rate is 0.040 (classic proteomics), 0.056 (transcriptomics) and 0.066 and 0.059 (metabolomics); the simulations gave 0.050, 0.056 and 0.044.
- **No RUV** is within the limits in all three, with a mean null rate of 0.04 to 0.05. The simulations gave 0.01 in the cells with a strong unwanted factor (C4) and 0.05 in the cell without one (C3: 0.050, 0.053, 0.048); the real figures sit at the second.

Things a reader should not miss, none of which crosses a limit:

- Metabolomics, recommended design, random splits: 13 of 100 splits have a call. That is above 0.10 and passes only because it is not above twice no-RUV's 12 of 100. Its mean null rate, 0.066, is the closest to the 0.075 limit of any passing row.
- On real data no RUV itself calls something in 4 to 12 of 100 null splits, against 0 to 1% in the simulations' global-null cell. The cause was not looked into.
- The per-split null rate ranges widely for every design (up to 0.2 to 0.5 in single splits). The limits are read against the mean, as revision 1.2 fixed.
- The batch-aligned splits of MTBLS404 gave the same picture as the random ones: tying the recorded batch to the false groups at a correlation of 1/3 did not raise any design's rate. So this data set did not put the designs under the condition where the simulations show today's design failing most (cell C5). The deposited table is the study's processed table; whether batch was already corrected in it was not verified at the source.

## Limits

- A null split measures false calls, not power (revision 1.1). Nothing here says the recommended design finds more.
- One data set per omic is one unit; 100 splits of it share samples.
- Kernel level, on the simulations' tree: the unmerged package kernels of the recommended design (`kit-126`, `kit-127`, `kit-129-k-dia`, `kit-129-trx`) are not what ran.
- Preparation is minimal and declared (no between-sample normalisation beyond what the deposited tables carry; the GUI route's filters and normalisation are not run).
- Metabolomics has 120 features, so a split's null rate moves in steps of 0.008.

## Not checked on real data

| Workflow | Data set accepted in revision 1.3 | State |
|---|---|---|
| DIA_LIMPA | PXD024043 | prepared for the runner, not run |
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
| DIA_LIMPA | random, batch | 0 | | | |

Laptop, one R process at a time, `nice -n 19`, one attempt each, 2026-11-18 10:26 to 11:18 (`runs/confirm.log`). Sidecars: `runs/<omic>/meta.txt` (tree commit, libraries, versions, command, seed rule, input sha256). The lifecycle test before them: `lifecycle/` (3 splits per omic and split type; the deliberate failure recorded as a crash and the run completed, in all four omics).

## To launch DIA_LIMPA on a hosted runner

From the `kit-125` work tree (the ticket pushes nothing):

```
git push origin kit-125:run-k125
```

This starts `.github/workflows/k125-null-split.yml` (one job; about 3 hours if the runner is as fast as the loaded laptop; the script's limit is 300 minutes): the library as `k124-sim.yml` builds it (needs the repository's `DEPS_READ_TOKEN`), the package tree `0b6fcb6d`, the PRIDE archive downloaded by the runner (658 MB, CC0) and checked against `data-manifest.tsv`, the lifecycle splits printed beside the laptop's, then splits 1 to 100 of both types. The artifact `k125-dia-limpa-1-dia-limpa` holds `dia_limpa/runs.tsv`; put that folder under `runs/dia_limpa/` and rerun `k125_aggregate.R` for its table. The workflow file is modelled on `k124-sim.yml` and has not been executed (no push from this ticket); its first run is also its test. One attempt: a failed job is recorded, not rerun. Runner figures sit beside the laptop's only after the lifecycle comparison is read under the CI tolerance rule.
