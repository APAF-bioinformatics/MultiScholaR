# KIT-151 cap evidence follow-up: results

Declaration `cap-followup-declaration.md` (12294f1a, pushed alone before the counted run). Lifecycle run 44653902387
(f8f4761b); counted run 44654405980 (a6f9ee12): 40 jobs, 3,060 of 3,060 rows, 0 failed cells, scored once by
`k151-followup-score.py` (sha256 20f4abbd..., unchanged; `k151-followup.R` 631959aa..., unchanged). No correction
pass. Hosted time 14,832 job-seconds on 4-core runners: about 16.5 core-hours (ceiling 30). The first run's result
(`cap-evidence-results.md`) stands as read.

Files: `followup-rows.tsv` (per scenario x seed x arm x mode x planted factor), `followup-summary.tsv`,
`followup-composition.tsv` and `followup-composition-summary.tsv` (kept and dropped features by what they carry),
`followup-jobs.tsv` (CPU and versions per job), `followup-score.json` / `followup-score.md` (readings),
`followup-source-excerpts.txt` (mofapy2 convergence lines), `followup-lifecycle-rows.tsv`.

## Recovery: seeds of 10 with |r| >= 0.9, F1 / F2 / F3

N = 20, weak F2 (the first run's data):

| arm | fast | medium | slow | slow iterations (median; seeds at maxiter) |
|---|---|---|---|---|
| none | 2/0/1 | 2/0/2 | 7/0/6 | 411; 0 |
| scale_views | 2/0/1 | 3/0/2 | 8/0/6 | 544; 0 |
| cap250 | 10/1/10 | 10/1/10 | 10/1/10 | 944; 5 |
| cap500 | 10/2/10 | 10/2/10 | 10/1/10 | 711; 2 |
| cap1000 | 10/2/9 | 10/1/10 | 10/1/10 | 731; 1 |
| cap2000 | 10/0/6 | 10/0/9 | 10/0/9 | 894; 1 |
| cap5000 | 6/0/3 | 8/0/8 | 8/0/8 | 1001; 8 |
| weight_views | 4/0/0 | 5/0/0 | 5/0/0 | 1001; 10 |
| zscore (workbook) | 2/0/2 | 5/0/4 | 8/0/8 | 489; 0 |

N = 20, strong F2:

| arm | fast | medium | slow | slow iterations (median; seeds at maxiter) |
|---|---|---|---|---|
| none | 2/0/2 | 2/1/2 | 6/1/6 | 391; 0 |
| scale_views | 2/0/1 | 3/1/2 | 7/1/6 | 549; 0 |
| cap250 | 10/9/10 | 10/9/10 | 10/10/10 | 966; 5 |
| cap500 | 10/7/10 | 10/7/10 | 10/7/10 | 846; 2 |
| cap1000 | 10/6/9 | 10/7/10 | 10/9/10 | 834; 2 |
| cap2000 | 10/6/6 | 10/6/9 | 10/7/9 | 939; 3 |
| cap5000 | 6/3/3 | 8/6/7 | 8/6/8 | 1001; 8 |
| weight_views | 3/7/0 | 2/7/0 | 2/6/0 | 1001; 9 |
| zscore | 2/0/2 | 5/2/3 | 8/2/8 | 496; 0 |

N = 60 (fast and slow only, as declared):

| arm | weak fast | weak slow | strong fast | strong slow |
|---|---|---|---|---|
| none | 8/0/4 | 10/7/10 | 9/5/5 | 10/9/10 |
| scale_views | 8/0/5 | 10/9/10 | 9/5/5 | 10/9/10 |
| cap250 | 10/10/10 | 10/10/10 | 10/10/10 | 10/10/10 |
| cap500 | 10/10/10 | 10/10/10 | 10/10/10 | 10/10/10 |
| cap1000 | 10/10/10 | 10/10/10 | 10/10/10 | 10/10/10 |
| cap2000 | 10/8/10 | 10/10/10 | 10/10/10 | 10/10/10 |
| cap5000 | 10/7/10 | 10/9/10 | 10/10/10 | 10/10/10 |
| weight_views | 10/8/8 | 10/9/8 | 10/10/8 | 10/10/8 |
| zscore | 7/0/4 | 10/9/10 | 7/3/4 | 10/10/10 |

Median iterations under fast: uncapped 21 to 34, capped 36 to 114, weight_views 106 to 219; under slow most arms run
400 to 1,000 and many capped / weight_views runs hit maxiter (1,000).

## Pre-declared readings

- (a) Early stopping: **No (not only).** On N = 20 weak under slow, none reaches F1 7 / F3 6 and scale_views
  F1 8 / F3 6, below the 9-of-10 mark. Slow training lifts the uncapped arms a lot (from 2/1 under fast) but not to the
  capped arms' level (10/10). At N = 60 slow, the uncapped arms do reach it (10/10). So at N = 20 the first run's
  poor uncapped result was partly early stopping and partly the imbalance.
- (b) Convergence default: **slow.** fast fails 10 comparisons (for example N = 20 strong cap1000 F2: 6 under fast vs
  9 under slow; cap2000 F3 6 vs 9). medium fails one (N = 20 strong cap1000 F2: 7 vs 9).
- (c) Cap: **1,000** (under slow). Per scenario, caps meeting F2 >= 9 and F3 >= 9: N = 20 weak none (F2 at most 1 of
  10 for any cap); N = 20 strong 250 and 1,000 (500: F2 7; 2,000: 7; 5,000: 6); N = 60 weak and strong all five.
  Largest cap meeting it in every discriminating scenario (N = 20 strong, N = 60 weak, N = 60 strong): 1,000.
  Because N = 20 weak does not discriminate while the others do, the pre-declared conclusion applies: **F2 at the weak
  strength is not detectable at N = 20 in this design by any treatment**, a limitation of the claim.
- (d) weight_views: **does not meet** the marks in any discriminating scenario under slow (N = 20 strong F2 6 / F3 0;
  N = 60 weak 9 / 8; N = 60 strong 10 / 8). It loses F3 (the large-only factor) at N = 20 in every mode and in 2 of 10
  seeds at N = 60. By the rule it is not offered.

## What the cap keeps (median over seeds; kept / total in that data set)

N = 20, weak data (D3, D4 are below every cap and are kept whole):

| cap | D1 F1-loaded | D1 F3-loaded | D1 unloaded | D2 F1-loaded | D2 unloaded |
|---|---|---|---|---|---|
| none | 2000/2000 | 2000/2000 | 16195/16195 | 200/200 | 1800/1800 |
| 250 | 138 | 140 | 22 | 93 | 158 |
| 500 | 256 | 250 | 69 | 129 | 372 |
| 1,000 | 430 | 440 | 232 | 168 | 832 |
| 2,000 | 710 | 720 | 703 | 200 | 1800 |
| 5,000 | 1227 | 1227 | 2718 | 200 | 1800 |

At N = 60 the enrichment is stronger (cap 1,000 keeps 564 F1, 559 F3 and 29 unloaded in D1). The cap does two things
at once: it shrinks D1 and D2 toward the small data sets (balance), and it keeps mostly loaded features (signal
enrichment: in D1 the loaded share rises from 19% to about 77% at cap 1,000 and N = 20). This design cannot separate
the two; F2 lives in D3 and D4, which no cap touches, so its gain comes from what the cap removes elsewhere.

## Inflated-noise sensitivity (caution)

10% of D1's features (2,000, all unloaded) given noise sd 3.0 instead of 1.5, so they are the most variable features.

| N = 20 | weak fast | weak medium | weak slow | strong fast | strong medium | strong slow |
|---|---|---|---|---|---|---|
| cap1000 (clean) | 10/2/9 | 10/1/10 | 10/1/10 | 10/6/9 | 10/7/10 | 10/9/10 |
| noisy cap1000 | 10/1/0 | 10/1/0 | 10/1/0 | 10/10/0 | 10/10/0 | 10/10/0 |
| cap2000 (clean) | 10/0/6 | 10/0/9 | 10/0/9 | 10/6/6 | 10/6/9 | 10/7/9 |
| noisy cap2000 | 10/0/9 | 10/0/10 | 10/0/10 | 10/5/10 | 10/6/10 | 10/7/10 |

With cap 1,000 the D1 selection is about 1,000 of the 2,000 inflated-noise features and almost no loaded feature
(median 0.5 F1-loaded, 0.5 F3-loaded kept): the large-only factor F3 is lost in every seed and mode. With cap 2,000,
about 1,830 inflated and about 100 F3-loaded features are kept, and F3 is still found (9 to 10 of 10). A variance cap
helps here only because the most variable features carry the signal; when unrelated features are more variable and
outnumber the cap, a factor that lives only in the large data set can disappear entirely.

## Instrumentation check (not a gate): reproduction of the first run

The fast-mode weak cells repeat the first run on the same data. 281 of 420 |r| values are identical, 21 differ by
more than 1e-6. 19 of those 21 (all differences above 1e-3) come from one job, N = 20 seed 151007, which ran on an
AMD EPYC 7763 64-Core Processor (`followup-jobs.tsv`, k151f-n20weak-7; the first run's job k151-main-7 ran on an
Intel Xeon 6973P-C, `jobs.tsv`); there iteration counts differ too (weight_views 166 in the first run, 106 here) and
two cells cross 0.9 (weight_views F1 0.854 in the first run -> 0.925 here, cap5000 F1 0.900 -> 0.903).
(Corrected 2026-11-20 on branch run-mofa-evidence-k151b: the earlier text had the two CPU models reversed; audit
SCI-05.) The remaining
differences are rounding-level on other CPU pairs. So double precision did not give bit-identical models across CPU
models here; this bears on INT-6 ("double: the same model on every computer") and is passed to the user, not acted on.

## Not as declared

Nothing. All cells ran; the reduced cross (N = 60 without medium) was fixed in the declaration. One note: MOFA2's
iteration count for a run that reached maxiter reads 1001 (the counter includes the start), and stopped-by uses
iterations >= 1000.

## Limits

One generator (Gaussian, linear, equal effects except the strong F2, one scale); N = 20 and 60; one training seed;
N = 60 without medium; one inflated-noise pattern; recovery depends on the CPU model at the margin (above).
