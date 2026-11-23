# KIT-183 lifecycle runs L1 (MOFA2 evidence lane, branch run-mofa-evidence-k183)

| run | package commit | what was run | observed |
|---|---|---|---|
| 45148992728 | 11545cb3 | 5 replicates per cell (n = 12, 20, 6; 12 workers, at most 4 at a time, one thread each), F3, the memory case, the declared failing replicate | every one of the 15 replicates `did_not_complete` ("$ operator is invalid for atomic vectors", no call stack); F3 identical; memory case 9.5 s, peak RSS 592 -> 1,176 MiB; declared failure (n = 7) recorded with its message; setup 633 s, simulation 980 s |
| 45151737782 | 51ca10b6 (the replicate records its call stack) | the same | the same 15 failures, located: `stats::aggregate(cbind(..., called = as.integer(planted$bh), ...) ~ ..., data = planted)`, where the data's column `planted` (the planted direction) masks the data frame `planted`; the failure is in the last aggregation, after training and every test, comparator and set of every data set and factor had run; F3 identical; memory case 10.3 s, peak RSS 592 -> 1,176 MiB; setup 649 s, simulation 986 s |

## Time per replicate (total_seconds, one thread, four workers at a time on a 4-core runner, Intel Xeon Platinum 8573C)

| cell | run 45148992728 (r1 to r5) | run 45151737782 (r1 to r5) |
|---|---|---|
| n = 12 | 63.6, 63.8, 64.0, 64.5, 53.4 | 59.9, 60.2, 59.1, 62.8, 50.4 |
| n = 20 | 65.4, 69.9, 65.9, 66.5, 42.6 | 61.1, 64.2, 60.8, 63.4, 38.5 |
| n = 6 | 65.6, 56.1, 58.0, 39.5, 32.1 | 61.6, 55.1, 54.9, 39.5, 32.3 |

The replicates that ran while four workers were busy (as in S) took 55 to 70 s; the last replicate of the workers that
ran two (fewer workers busy) 32 to 53 s. Training and analysis were not timed separately (the split is written only for
a completed replicate). Inferred: every replicate's time includes all its work, since the failure came after the last
test.

## Stop rule (ticket 11)

"If L1 measures more than 50 s per replicate on one core ... stop and return the numbers to the owner with a re-split;
the replicate count is not reduced to fit." Measured 55 to 70 s under S's load: **the rule fires; S is not started.**

Re-split for the owner (arithmetic, at 65 s per replicate and four workers per job): 300 replicates per job need
300 x 65 / 4 / 60 + 10 = 91 minutes, above the 75-minute timeout. 12 jobs of 200 replicates need 64 minutes each
(timeout 75 holds): S expected 12 x 64 + 13 (determinism) = 781 job-minutes = 52.1 core-hours (planned 26.08),
ceiling 12 x 75 + 20 = 920 job-minutes = 61.3 core-hours (planned 41.33); with S' at the same ceiling the ticket's
ceiling becomes 222.67 + 2 x (61.3 - 41.33) = 262.6 core-hours, above the 223 the owner accepted (VAL-53 (7)).

## Other L1 items

- Pool counts and sets formed (6.2): **not measured** (written only for a completed replicate).
- Memory case (8.7): 2,000 sets on 20,000 features x 500 samples in 9.5 and 10.3 s (ceiling 60 s); process peak RSS
  1,176 MiB (the ticket's 1 GB peak not met as worded, VAL-85 (5)); the high-water mark before the case was 592 MiB,
  so the case added at least 584 MiB.
- MOFA2's sample-count warnings: not recorded (written only for a completed replicate).
- The defect is corrected in `simulate-k183.R` after run 45151737782 (column names only inside the formulas); not run.

## L1 check after the aggregation fix (run 45219821547, package 874f6c45; VAL-54 instrument rule, VAL-64)

Observed: 15 of 15 replicates completed. Time per replicate under four workers: n = 12 57.4 to 60.5 s (training 20 to
22 s, analysis 35 to 38 s), n = 20 59.2 to 61.6 s, n = 6 53.4 to 60.0 s; the last replicates of two-replicate workers
31 to 49 s. Within the re-split's 65 s (VAL-64). F3 identical; the memory case as before; the declared failure recorded.
MOFA2's messages at n = 6 and 12 (n = 12 counts 11 biological samples in every view, the transcriptomics view lacking
one): "only 6 (11) biological samples are measured in every view (fewer than 16); factor models need more samples",
"The total number of samples is very small for learning 5 factors. Try to reduce the number of factors to obtain
meaningful results. It should not exceed ~3.", at n = 6 also "Some group(s) have less than 10 samples, MOFA will have
little power to learn meaningful factors for these group(s)".

Pools (means over the 5 replicates): proteomics complete features 1,217 (n = 12), 1,118 (n = 20), 1,344 (n = 6), with
494, 633, 200 ranked features incomplete; metabolomics complete 164, 151, 180, identifiers in U_cam 144, 133, 155;
transcriptomics 4,000 (no censoring). Classes with fewer than 90 % of their declared sets formed (16 class-size-cell
cells): proteomics C2 size 30 (n = 20: 0 %), C4 size 100 (n = 12: 60 %, n = 20: 0 %), P-up and P-down size 100
(n = 20: 80 %); metabolomics P-up and P-down size 10 (40 to 80 %) and 20 (0 %). Redefined in the declaration's dated
section before S: proteomics C2 sizes 10, 15; proteomics C4 and P sizes 10, 30; metabolomics P size 5.

Not reported (ticket 6.2 asks the lifecycle for them): the within-block residual correlation and the variance-trend
fit. They cannot be added without a run, and their absence is not a reason for one (AUDIT-KIT-183-impl I07). Run
times in the tables above are times since each job's start; the GitHub job times are in `runs.tsv` (L1 16.5, L1b 16.6
job-minutes).
