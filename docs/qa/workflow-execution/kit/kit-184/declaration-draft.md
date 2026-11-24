# KIT-184 simulation: declaration, DRAFT (third version)

**This is a draft. It is not the declaration, it fixes nothing, and nothing may be counted against it.** The counted
run S has not started and does not start from this file. This version applies the declaration audit
(`tickets/kit/audits/AUDIT-KIT-184-declaration-20261124.md`, findings K184-D01 to D12), the owner's decisions O1 and
O2 (methods ledger VAL-150, 2026-11-24) and the audit's re-check 1 (K184-R01 to R05, with its "order before the
freeze"). The declaration itself is then `declaration.md`, committed alone, before S, after the auditor has confirmed
the blind lifecycle extension L1b. Earlier versions of this draft: e9738ade (first), fc35c308 (second).

Ticket `tickets/kit/KIT-184.md` revision 1.7, sections 6.2 to 6.5 and 11; ledger VAL-11, VAL-53, VAL-58, VAL-62, VAL-94
(KIT-183's outcome O1), VAL-97 (3), VAL-150.

Where a paragraph below is the audit register's text it is used as the register gives it. Section 6 lists the points
on which the register's first texts and the code disagreed, how re-check 1 settled them, and two small additions I
made to L1b's fixed output list.

## 1. What I had seen of L1 when I proposed the changes

The changes of section 4 alter what the simulation measures, and they were proposed after L1's results had been
seen. This is everything of L1 (run 45333818847, 5 replicates of each of cells A, B and D) that I looked at. The list
was taken from the session's record of my commands, which is not in the repository; the counts section 4 rests on are
now on the branch (`l1/formed.tsv`, `l1/clean-features.tsv`) so that they can be checked without it.

**Rejection rates and other results I saw** (all from 5 replicates per cell; all in `l1/` on the branch):

- `l1/results-summary.md`, every line (cut at 900 characters): how many of K-E1's 33 checks passed in cells A and B
  (17 and 18) with the largest upper limit; K-E3's counts; K-E4's rate and upper limit per data set in cells A and B
  (in cell B all three data sets fail it: 0.080, 0.107, 0.200); K-E5's pass counts (19 and 21 of 33); the outcome rows.
- `l1/error-rates.tsv`: class C5 only, test `camera`, arms 1, 2, 3 and 3o, every cell, on F_main, F_ind and F_ret,
  pooled over set sizes (the table in `lifecycle.md`, and the same for F_ret). For proteomics C5 under arm 3 that is
  0.000, 0.107, 0.020 on F_main and 0.093, 0.120, 0.147 on F_ind (cells A, B, D). No rate by set size exists in L1's
  outputs: the simulation summed the sizes before writing.
- `l1/feature-rates.tsv`: cell A on F_ind and F_ret, every arm; arm 3 on F_main in every cell. Six rows of one
  worker's `feature.tsv` with their counts of significant features (cell A, replicate 1, F_main).
- `l1/power.tsv`: arm 3, `camera` and `camera_floor`, sets of 10, F_main, cells A and B (what [PF] and [PC] read).
- `l1/numbers.tsv`, in full, twice: the nine text numbers as 5 replicates would give them, with where each was found.
- `l1/caught.tsv`, `l1/identity.tsv`, `l1/exactness.tsv`, `l1/completion.tsv` in full; the first rows of `l1/floor.tsv`;
  the first rows of one worker's `counts.tsv` (two values of `factor_r2`); seconds per replicate.

From this I wrote, in `lifecycle.md` and in my report to the coordinator, an indication that K-E4 may fail in cell B
(arm 3 calls C5 sets on the main factor in 8 to 20 %, the oracle arm in 0 to 1 %).

**What I did not look at:** any rate of C4, C4b, C1-out or C3; the rows of `ke1-checks.tsv`, `ke4-checks.tsv` and
`ke5-checks.tsv` (only the summary's lines); `paired-differences.tsv`; `operating-characteristic.tsv`; the rest of
`power.tsv`; any per-replicate rate (`errors.tsv` of a worker: its line count only).

**Mechanics used for section 4** (counts only): `l1/pools.tsv` and the workers' `formed.tsv` (now `l1/formed.tsv`:
sets wanted and formed per replicate, model, data set, factor, class and size); the `features` column of the workers'
`feature.tsv` (now `l1/clean-features.tsv`: how many clean complete features entered the feature-level rate, per
replicate, data set and factor); the generator's constants in `simulate-k184.R`. The script did not write pool sizes:
where a pool size is given it is arithmetic from the constants or from those counts, or the auditor's
re-implementation of the pool logic (register section 1; not recomputed by me), and says which.

**Lifecycle replicates are excluded from S.** L1 ran replicates 1 to 5 of each cell, and the generator's seed depends
only on the cell and the replicate number. S therefore runs replicates 6 to 805 of cells A and B and 6 to 405 of cell D
(800, 800 and 400 replicates, none of them seen), and the determinism job replicates 6 to 15 of cell A. No L1 output
enters any count, rate or text number. Checked: over S's ranges no generator, F_ind or set-drawing seed is used twice.
The permutation seed 20261185 + r (KIT-183's comparator rule, reported arms only) recurs across cells at the same r
and, for r = 661 to 805, equals the cell A generator seed of replicate r - 655; each permutation test stays valid, its
permutations being independent of its own replicate's data. (The check is the auditor's, register Q2; the seed
formulas are the script's.)

## 2. Departure from the programme's blind-lifecycle rule

**Departure from the programme's blind-lifecycle rule.** The programme scores lifecycle runs blind (VAL-11). L1 was
not: the scorer committed before L1 ran on its five replicates per cell against the planted truth, and its tables are
committed in `l1/`. The ticket's L1 row asked for C5's rate under arms 1, 2, 3 and 3o (plan audit A02); the other
tables went beyond that row. What was read is listed in section 1. Consequences: no L1 replicate enters any count, rate
or text number; no cell, arm, test, limit, criterion or replicate count was changed after L1; the one change to a
deciding check is the owner's decision, ledger VAL-150; the lifecycle extension L1b is blind (section 5). L1's tables
stay committed, labelled as lifecycle output.

The owner has acknowledged this disclosure and decided to proceed without a re-run (O2, VAL-150).

For `results.md`, first section, when it is written: "Lifecycle disclosure: L1 (45333818847) was scored on real truth
at five replicates per cell before this declaration; what was seen and what was changed afterwards are in
`declaration.md` sections 1, 2 and 4."

## 3. The declaration as I would freeze it

### Vehicle, cells and jobs

- MOFA2 evidence lane `.github/workflows/k184-evidence.yml` on branch `run-mofa-evidence-k184` (the lane's
  `.github/mofa-smoke/k184.R` is a byte copy of `simulate-k184.R`, checked by the job): R 4.6.1, the KIT-134 snapshot
  library, MultiScholaR at the package commit named in the lane's `k184-params.json`, MOFA2 1.22.0, ruv 0.9.7.1,
  RUVIIIC 1.0.19, limma 3.68.4 (as L1 found them, `l1/versions.tsv`), one thread per worker
  (`OPENBLAS_NUM_THREADS = OMP_NUM_THREADS = MKL_NUM_THREADS = 1`), four workers at a time per job. Every replicate
  writes the versions of ruv, RUVIIIC, limma and MOFA2 it ran with; S's `versions.tsv` holds their unique rows.
- Scripts: `simulate-k184.R` and `score-k184.R` at that commit, their SHA-256 written here at the freeze, after L1b.
  KIT-183's two files are not edited.

| Cell | Biological samples | rho | k removed | Replicates | Role |
|---|---|---|---|---|---|
| A | 12 (transcriptomics 11: the last sample is missing) | 0 | 2 | 800 (numbers 6 to 805) | deciding |
| B | 12 (transcriptomics 11) | 0.6 | 2 | 800 (6 to 805) | deciding |
| D | 12 (transcriptomics 11) | 0 | 4 | 400 (6 to 405) | reported only |

Cell C is not run: KIT-183's counted outcome is O1 (`camera` is the test of record at every size; VAL-94).

- Jobs (re-split after L1, accepted by the coordinator in the form 16 jobs of 125 replicates): jobs s01 to s08 each
  hold 100 replicates of cell A (6 to 105, 106 to 205, ..., 706 to 805) and then 25 of cell D (6 to 30, ..., 181 to
  205); jobs s09 to s16 each hold 100 replicates of cell B (same ranges) and then 25 of cell D (206 to 230, ..., 381 to
  405). Timeout 85 minutes each; no replicate starts in the last 6 minutes (recorded `not_started`). Determinism job:
  cell A replicates 6 to 15 again, timeout 25 minutes.
- In every job the lane starts its four cell A or B processes (25 replicates each) at once and a cell D process (6
  or 7 replicates) only when one of them has finished, at most four at a time, so a job that runs short loses
  replicates of cell D (reported only) first. At L1's largest times a worker's last replicate starts at about 62
  minutes against 68; on a runner about 10 % slower than L1's, the last cell D replicate of each worker would not start
  (up to 64 of 400; cell D 'not decided', R6). Each job records its runner image and processor model.
- A replicate in which MOFA2, the correction or a test fails is `did_not_complete`. A cell with more than 5 % of its
  declared replicates not completed (failed or not started) is "not decided".

### Generator (per replicate r of a cell)

One stream, `set.seed(20261840 + offset + r)` with offset 0 (A), 100,000 (B), 300,000 (D). Everything of KIT-183's
declaration (`kit-183/declaration.md`, "Generator") is unchanged: features (4,000; 1,800; 250 on 204 identifiers),
abundances, noise, Z1 on 10 % of the loading units and Z2 on another 5 % of them (none in metabolomics),
value-dependent missing values in proteomics and metabolomics, the metabolomics identifiers carrying two or three
features. Added:

- Per data set two unwanted factors over the biological samples, independent between data sets: U2 ~ N(0, 1);
  U1 = rho g + sqrt(1 - rho^2) e, g the group indicator scaled to mean 0 and variance 1 in the sample, e ~ N(0, 1) made
  orthogonal to g and scaled in the sample, so the sample correlation of U1 with the group is exactly rho.
- Block candidates: before any value is drawn, a random share of the features with zero loading on Z1 and Z2 on
  single-feature identifiers is set aside and carries no unwanted loading: transcriptomics 0.30, proteomics 0.50,
  metabolomics 0.55. (The ticket leaves these shares to be fixed from L1's counts so that every block forms in at least
  90 % of replicates. In L1 every block formed in 15 of 15 replicates in every data set, so the indicative shares stand.)
- Every other feature carries U1 with probability 0.4 and U2 with probability 0.4, independently; a carried loading has
  random sign and size U(1, 2) times the feature's noise standard deviation.
- Values y = mu + sigma (l1 Z1 + l2 Z2 + lu1 U1 + lu2 U2 + e). The mask is drawn on those values, before the block
  latents. One run per biological sample.
- Null blocks from the complete candidates, disjoint, in the order rho_block 0.1, 0.3, 0.6, then size, then block:
  transcriptomics 2 blocks per rho_block and size (sizes 10, 30, 100), proteomics 1 (10, 30, 100), metabolomics 1 (5, 10).
- "Clean": zero loading on Z1, Z2, U1 and U2. "Complete": no missing value.

### The correction, the inputs and the models

- Per data set RUV-III with the two groups as replicate sets, every feature a potential control and k factors:
  transcriptomics `ruv::RUVIII(Y, M, ctl = every gene, k)` on its 11 samples; proteomics and metabolomics
  `RUVIIIC::RUVIII_C_Varying(k, Y, M, toCorrect = every feature, potentialControls = every feature)`, then
  `restoreRuvFeatureLevels()`, rows without a finite corrected value dropped.
- Inputs written by `writeIntegrationInput(..., corrected = )` with `correction = list(state = "simulated", method,
  k, replicate_structure = "group", grouping_variable = "group")`: the removed factors are the writer's own.
- Two models per replicate, `runIntegrationModel()` at every shipped default except `training_seeds = 1`: one with
  `input_matrix = "ruv_corrected"` ("the corrected model"), one at the default input ("the default model").

### Factors and covariates

- Corrected model: **F_main** (largest absolute correlation of the oriented scores with Z1). **F_ind**: N(0, 1) per
  biological sample, seed 20761840 + r (A), 20861840 + r (B), 21061840 + r (D); no model's factor. **F_ret** per data
  set: F_ind projected onto the space that data set's corrected values span (right singular vectors with singular value
  above 1e-8 of the largest, of the row-centred corrected matrix of the features complete in that data set's samples).
- Default model: its F_main, and **F_unw** (largest absolute correlation with any planted unwanted factor).
- Per data set and factor: `integrationEnrichmentRank()` (Pearson, 5 samples), `integrationEnrichmentMappedList()`,
  `integrationEnrichmentTestValues()`; sets drawn with seed 20261888 + (generator seed modulo 1,000,000) + 7919 x the
  running number of the model, data set and factor.

### Arms and tests (the same sets under every arm of a model, data set and factor)

| Arm | Model | Values | Design | Tests |
|---|---|---|---|---|
| 3 | corrected | before RUV | intercept, factor, removed factors | `camera`, `camera_pr`, `camera_floor` |
| 2 | corrected | before RUV | intercept, factor (`adjust_ruv_factors = FALSE`) | `camera`, `camera_pr`, `camera_floor`, `gsea`, `permutation` |
| 1 | corrected | corrected | intercept, factor (script only) | `camera`, `camera_pr` |
| 3o | corrected | before RUV | intercept, factor, the planted U1 and U2 (script only) | `camera`, `camera_pr`, `camera_floor` |
| D0 | default | before RUV | intercept, factor (setting absent; `design_covariates` must be `none`) | `camera`, `camera_pr`, `gsea`, `permutation` |
| D1 | default | before RUV | intercept, factor, removed factors (`adjust_ruv_factors = TRUE`) | `camera`, `camera_pr`, `camera_floor` |

Arms 1, 2, 3, 3o on F_main, F_ind and F_ret; D0 and D1 on the default model's F_main and F_unw. Metabolomics on a
learned factor goes through `runIntegrationEnrichment()` with a synthetic KEGG table (the whole product path;
`camera_pr` through the kernel); everything else through `integrationEnrichmentPathwayTest()`. Sizes 5 to 500, BH
within the call (one call per model, data set, factor, arm and test, over every class's sets). `gsea` seed 20261121;
`permutation` 1,000 reorderings, seed 20261185 + r. `camera_floor` and `permutation` are reported only (VAL-97 (3)).
Under every arm that runs both, p_floor at least p_camera is checked set by set.

### Classes (members are identifiers of the test's universe; drawn per replicate, model, data set and factor)

What the simulation draws and tests is the ticket's, unchanged, as the scripts committed before L1 draw it:

| Class | Members | Sets per replicate | Sizes drawn: transcriptomics; proteomics; metabolomics |
|---|---|---|---|
| C0 | uniform from the universe | 10 per size | 10, 30, 100; 10, 30, 100; 5, 10, 20 |
| C1 | each formed null block, by rho_block | 2; 1; 1 per rho_block and size | 10, 30, 100; 10, 30, 100; 5, 10 |
| C2 | the lowest-abundance fifth (metabolomics: third) of the clean complete features in no block | 5 per size | 10, 30; 10, 30; 5 |
| C5 | complete features that are not candidates, zero loading on Z1 and Z2, positive loading on U1 | 5 per size | 10, 30, 100; 10, 30, 100; 5 |
| P-up, P-down | half with a positive (negative) Z1 loading, half clean complete in no block | 5 per size and direction | 10, 30, 100; 10, 30; 5 |
| C4 | clean complete features in no block that the model selected | 5 per size | 10, 30, 100; 10, 30; 5, 10, 20 |
| C1-out, C4b | block members the model did not select (at least 5; not metabolomics); those it selected (at least 5) | per block | as the blocks |
| C3 | metabolomics identifiers carrying two or three zero-loading features | 5 per size | 5, 10 |

A set of a size forms in a replicate when its pool holds at least that many features; sets of one class and size are
drawn independently from the same pool. The simulation writes its counts per set size; every redefinition after L1 is
made in the scorer (`score-k184.R`, `CRITERIA184`), which sums the sizes before every estimand:

| What the scorer reads differently from the sets drawn | Where |
|---|---|
| Proteomics C5 in two compositions: sizes 10 and 30; and every size formed (10, 30 and, where it formed, 100) | "Proteomics C5" under the criteria below; section 4.1 |
| Reported C4 rate, transcriptomics: sizes 10 and 30 under the corrected model; not reported under the default model | section 4.2 |
| Reported C4 rate, metabolomics: sizes 5 and 10 | section 4.3 |
| Feature-level rate on F_main in metabolomics: not defined | section 4.4 |

Expected formation in S, from the generator's arithmetic: every block, C0, C1, C2 and P size in every replicate;
proteomics C5 sizes 10 and 30 in every replicate; metabolomics C5 (size 5) in about 93 % of replicates; the last
metabolomics block in about 99 %. T is the mean over the replicates with a set of the class, and the number of those
replicates is reported with every T; for metabolomics C5 the rate is conditional on at least five complete U1-positive
features. A class that forms in fewer than 90 % of S's replicates is not redefined after S; it is reported with its
count.

(The 93 % and 99 % are the auditor's re-implementation of the pool logic, 20,000 replicates, register D07 and section
6; L1 formed metabolomics C5 in 15 of 15. The consequence: the metabolomics C5 check of K-E1, deciding, is expected to
rest on about 750 of the 800 replicates of a cell, those in which at least five such features are complete.)

### Estimands (as ticket 6.3)

- Error rate: s_r = the share of the replicate's sets of the class (all sizes together) with unadjusted p at most
  0.05; T = the mean of s_r over the replicates with a set; interval T plus or minus t(R - 1) sd(s_r) / sqrt(R).
- Paired differences T1 minus T3 and T2 minus T3 per class, data set and factor, on F_ind, F_ret and F_main.
- Feature-level rate per arm (1, 2, 3, 3o): the share of clean complete features with moderated-t p at most 0.05 for
  the factor, on F_ind, F_ret and F_main; on F_main over the clean complete features the model did not select (not
  defined for metabolomics, section 4.4).
- Power: the share of P-up and P-down sets with BH-adjusted p at most 0.05 in the expected direction on F_main, per
  arm, test, size, direction and data set; direction accuracy with the scorer's own orientation.
- Exactness per data set: the status of the written factors, the largest s_(k+1)/s_1, the largest relative residual.
  The identity check per replicate, data set and factor (arm 3's `correlation` against the same model on the corrected
  values; below 1e-8). What the removed factors caught (smallest canonical correlation with the planted factors).
- Determinism: replicates 6 to 15 of cell A run again in a separate job; twelve tables must be identical to the same
  replicates in S: `errors`, `power`, `formed`, `counts`, `feature`, `identity`, `floor`, `exact`, `caught`, `factors`,
  `orientation`, and `status` without its seconds columns. A difference is reported with its size and cause (plan
  section 4); a difference in a table that a criterion reads is a defect, handled as R3 within S' (275 job-minutes);
  one that needs more is a stop to the owner.

### Criteria, outcomes and the numbers in the texts (as ticket 6.4 and 6.5; nothing added, no threshold changed)

- **K-E1** (decides in A and B): arm 3, `camera`, upper 95 % limit of T at most 0.075 on F_main and F_ind for C0, C1
  (each rho_block) and C2, and on F_ind for C5, in every data set: 33 checks per cell.
- **K-E2** (`camera_pr`, C0 and C2): decides nowhere under KIT-183's O1; run and reported.
- **K-E3**: direction accuracy of arm 3 at least 0.99 in cells A and B, decided when at least 50 planted sets were called.
- **K-E4** (declared, does not hold the merge): arm 3, `camera`, the same limit on C5 on F_main, every data set, A and B.
- **K-E5**: arm 2, `camera`, the same limit on K-E1's 33 checks.

**Proteomics C5 (dated, after L1).** Size 100 formed in 10 of L1's 15 replicates (by the generator's arithmetic in
about 55 to 65 % of replicates), below the 90 % rule of ticket 6.2. The rule's redefinition (sizes 10 and 30) was
chosen after L1's pooled C5 rates had been seen (section 1). Every check, number and outcome that reads proteomics C5
therefore reads it in two compositions: (i) sizes 10 and 30; (ii) every size formed in the replicate (10, 30 and, where
it formed, 100), as the ticket declared and the scorer committed before L1 computes. A check (K-E1 and K-E5 on F_ind,
K-E4 on F_main) holds for proteomics C5 only when it holds in both; [R3], [R5] and [R6] take the larger T of the two;
[R5m] is arm 3's T in the composition and cell in which [R6] was found. Size 100 alone, over the replicates in which it
formed, is reported with that count; it is conditional on a pool of at least 100 complete features, which is not
independent of the data, and decides nothing. (The owner's decision O1, methods ledger VAL-150, 2026-11-24.)

- Outcomes R1 to R8 exactly as ticket 6.5. An unfavourable valid result is reported, not re-run. R6 is read as: a
  cell "not decided", or the written factors not `exact` in more than 1 % of a data set's replicates, or the identity
  check above 1e-8 in any replicate, read from `completion.tsv`, `exactness.tsv` and `identity.tsv` (the scorer's
  summary has a line for it).
- [R3], [I3], [R5], [I5], [R6], [I6], [R5m], [PF], [PC] as ticket 6.4 and 5.5 define them, with the proteomics C5
  paragraph above, as the scorer's `numbers.tsv` gives them.
- Reported, deciding nothing: arms 1, 2, 3o, D0, D1 on every class; C4, C4b, C1-out, C3; the rates per set size
  (`error-rates-by-size.tsv`, each with the number of replicates in which the size formed); the paired differences; the
  feature-level rates; power; cell D; `camera_floor` and `permutation`; the exactness measures; the operating
  characteristic from the measured standard errors. `operating-characteristic.tsv` counts the proteomics C5 check
  twice (34 rows per cell), which can only understate the probability of passing.
- C4, C4b and C1-out are defined by the model's selection, so which features and which replicates enter them depends
  on how much unwanted variance the correction left; their rates are conditional on that and are not calibration
  evidence (KIT-183 measured C4 at about 10 % on an independent covariate, VAL-94).

### Attempts

One counted run S. One correction S' only for a failure of S's instrumentation, limited to re-running the affected
jobs within what the ceiling leaves (305 job-minutes less the score job's timeout: 275); if a whole repeat were needed,
stop, to the owner (ticket 11, revisions 1.5 and 1.7). An unfavourable valid result is recorded and its outcome
applied.

The scorer runs once, in the lane's `score` job after every evidence job, with a timeout of 30 minutes (three times
L1b's measured 119 seconds is 6 minutes, below the floor of 30), on the artifacts, which are kept 90 days. Its minutes
count in row S: expected about 3, ceiling its timeout, taken from S', which becomes 305 minus that timeout, 275
job-minutes. A score
job that fails or times out has produced no result: the scorer's mechanics are fixed without changing any estimand,
check, limit or composition, and the same artifacts are scored again (VAL-54, at most three). S is never re-run for it.

### Dated section (after L1, before S), as it would read

- Lifecycle: run 45333818847 (15 of 15 replicates completed; parity with KIT-183's cell identical; 110 to 133 s per
  replicate under four workers). The stop rule of ticket 11 fired at 70 s; the coordinator accepted the re-split above
  and kept the ceiling (ticket 11, revision 1.5). The first push (45333717225) ran no job. L1 was not scored blind
  (section 2). The lifecycle extension L1b (run 45348572248; every item of the list of section 5 held, `lifecycle.md`;
  to be confirmed by the auditor before the freeze) was blind (section 5).
- Frozen with: `simulate-k184.R` SHA-256 58333b1ef1b479a65f30089dd0307aad5ebab2495c1fb060621a8738a873a75b and
  `score-k184.R` 282b844e1c3427e8c99f296aa486c313957776bda8067d782a9a4511b43fa8b2 (package commit 9e525fd9, the commit
  L1b ran); lane commit ef6ce13d2b965a8cbaf7c0280ba941498c3a4273, `k184-evidence.yml` SHA-256
  ae0dc5a52c2735d5b6092fca60455c4f1e9093a7de08b3a7b86c98f04e03eea7 (the workflow L1b ran, 2b680757..., with the score
  job's timeout set to 30 minutes and nothing else changed); the score job's timeout 30 minutes.
- After L1's formation counts (section 4): proteomics C5 read in two compositions, both of which must hold (VAL-150);
  the reported C4 rate without the sizes that did not form; the feature-level rate on F_main not defined for
  metabolomics. All in the scorer; the sets drawn and tested are the ticket's. K-E1 keeps 33 checks per cell and its
  limit.
- S runs replicates 6 onward; L1's and L1b's replicates enter nothing.

## 4. Changes forced by L1

The ticket's rule (6.2, last sentence; KIT-183's rule): "a class with fewer than 90 % of its declared sets formed is
redefined in a dated section of the declaration before the counted run."

The rule is applied as the scorer committed before L1 counts it (`score-k184.R` at 21053957, lines 394-398): per
cell, model, data set, factor, class and size, a size being redefined when it is below 90 % in any cell; KIT-183's
dated section applied it the same way. The ticket's sentence names the class; read per class, the transcriptomics C4
size 100 under the corrected model and the metabolomics C4 size 20 would stay (each 70 of 75 sets in cell B) and
proteomics C5 would still be redefined. With five replicates per cell the rule flags a size that forms in 93 % of
replicates with probability 0.66 (the generator's arithmetic puts metabolomics C4 at size 20 at about 93 %); no
deciding check depends on the reading.

The coordinator named three points; I found the fourth (4.3) while writing the first version. No size, class,
criterion, margin, cell or arm is added and no threshold changed. Every change is made in the scorer: the generator and
the set drawing stay as the scripts committed before L1 have them (register K184-D03), and the only change to the
simulation script is that its error table now has one row per set size.

### 4.1 Proteomics C5 at size 100: decided by the owner (O1, VAL-150)

- **Declared:** C5 in proteomics at sizes 10, 30 and 100, 5 sets per size. The estimand is the per-replicate share
  over whatever sets formed: 15 where size 100 formed, 10 where it did not. The ticket's own indicative pool is "about
  100" features (6.2), written before any run.
- **L1's mechanics:** size 100 formed in 10 of 15 replicates (4, 3, 3 of 5 in cells A, B, D), identically under both
  models, since the pool does not depend on the model. Sizes 10 and 30 formed in 15 of 15.
- **Arithmetic** (generator constants): 1,800 features, 180 with a Z1 loading, 90 with a Z2 loading, 1,530 with
  neither; half of those (765) are candidates; of the other 765 a share 0.4 x 0.5 = 0.2 has a positive loading on U1,
  about 153; at the ticket's indicative completeness of 0.674 about 103 are complete (standard deviation about 9.5),
  which holds 100 features in about 65 % of replicates. The auditor's re-implementation, which draws the mask on the
  values as the code does, gives a pool of 101.1 (standard deviation 9.4) and 56 %: features loaded on U1 are complete
  a little less often. L1 found 67 %. Where size 100 forms, its five sets are five draws of 100 from about 100 to 115
  features, so they are nearly one set counted five times.
- **Could the choice have been steered?** I proposed dropping the size having seen this check's pooled rate above the
  limit in every cell (section 1), without any rate by size. The auditor's reading from the design alone (register
  section 1 (c)): on the main factor, wherever the correction leaves unwanted variation aligned with the factor, larger
  sets of features that share it are called more often, so dropping size 100 can only lower the proteomics C5 rate
  there, the favourable direction for K-E4, [R5] and [R6]; on the independent covariate, which decides, the direction
  cannot be known from the design. That is why the remedy was the owner's to choose.
- **The four options, and what each is in the scorer** (`CRITERIA184$c5_proteomics`, one line):

| Option | The scorer's line | K-E1 and K-E5 (proteomics C5 on F_ind) | K-E4 (proteomics C5 on F_main) | [R3], [R5], [R6], [R5m] |
|---|---|---|---|---|
| (a) sizes 10 and 30 only (the first draft) | `list(sizes_10_30 = 100L)` | one row, 10 sets per replicate | one row | read from that row |
| **(b) both compositions must hold (decided)** | `list(sizes_10_30 = 100L, declared = integer())` | two rows; the check holds when both hold; still 33 checks per cell | two rows; holds when both hold | the larger T of the two; [R5m] in [R6]'s composition |
| (c) the declared composition only | `list(declared = integer())` | one row, 15 or 10 sets per replicate | one row | read from that row; the 90 % rule set aside |
| (d) a lower proteomics candidate share so that size 100 forms | as (c), and the generator changed | as declared, on other data | as declared, on other data | as declared; needs another lifecycle run |

  For (d) the auditor's re-implementation gives every block formed / size 100 formed in 0.977 / 0.943 of replicates at
  a share of 0.43, 0.997 / 0.922 at 0.44 and 1.000 / 0.887 at 0.45: the window is about 0.43 to 0.44, not the 0.43 to
  0.45 the first draft gave. It would change the data under every proteomics check after every pooled L1 rate was seen.
- **Decided: (b).** It cannot be easier than what was fixed before L1, it meets the rule, and the generator is not
  touched. The text is "Proteomics C5 (dated, after L1)" in section 3.
- **What it changes:** one of K-E1's 33 checks per cell (and the same check of K-E5, and K-E4's proteomics check) is
  read twice and must hold both times. The limit and the number of checks are unchanged. Its worst-case Monte Carlo
  standard error is unchanged (0.0077 at 800 replicates and a true rate of 0.05); the measured one may be smaller or
  larger, since the five sets of 100 were near-copies that carried a third of a replicate's weight. A calibrated test
  fails K-E1 by chance slightly more often than with one composition (register section 5: at most the worst-case 0.10
  of one more check, in practice far less because the two compositions share their sets of 10 and 30).

### 4.2 Transcriptomics C4

- **Declared:** C4 ("as KIT-183 forms them", on the restated pool): clean complete features in no block that the model
  selected; transcriptomics sizes 10, 30, 100; 5 sets per size. Reported for every arm; decides nothing; no text number
  reads it.
- **L1's mechanics:** under the default model size 10 formed in 2 of 15 replicates and sizes 30 and 100 in none. Under
  the corrected model sizes 10 and 30 formed in 15 of 15 and size 100 in 14 of 15 (4 of 5 in cell B).
- **Arithmetic:** the model keeps 1,000 of the 4,000 transcriptomics features, those with the largest variance above
  the trend. About 1,877 features are clean (the 1,020 candidates and 0.36 of the other 2,380 zero-loading features; L1
  counted 1,841 to 1,926 per replicate, `l1/clean-features.tsv`, F_ind rows), of which about 1,037 are in no block.
  Default model (values before RUV): about 1,907 features carry at least one unwanted loading of 1 to 2 noise standard
  deviations, more than the 1,000 the model keeps, so clean features are almost never kept and the pool is nearly
  empty by construction. Corrected model, size 100: not from arithmetic alone, since it depends on what the model
  selects after the correction. From L1's counts (the F_ind row less the F_main row of `l1/clean-features.tsv`) the
  model kept 306 to 362 clean features in 14 replicates, of which about 0.55 (170 to 200) are in no block; in
  replicate 4 of cell B it kept 77 (about 43 in no block), so size 100 did not form there.
- **Options:** (a) per model: the reported rate over sizes 10 and 30 under the corrected model; no transcriptomics C4
  rate under the default model (arms D0 and D1). (b) One rule for both models: every size is flagged somewhere, so no
  transcriptomics C4 rate at all. (c) Another pool under the same name (zero loading on Z1 and Z2 whatever the unwanted
  loadings): a different class. (d) Under the per-class reading of the rule the corrected model would keep size 100.
- **Proposed: (a)**, in the scorer.
- **What it changes:** Nothing the simulation draws or tests changes: the sets are drawn as the scripts committed
  before L1 draw them (`simulate-k184.R` keeps C4 at its declared sizes), so every p-value, BH-adjusted p, power
  figure and K-E3 count is what the pre-L1 code gives. The scorer leaves the size out of the reported C4 rate and
  reports, for each size left out, how many replicates it formed in. (The first draft removed the sizes in the
  generator and said "no criterion, no check, no text number": that was false. Every set of a data set, factor, arm
  and test is in one Benjamini-Hochberg family, and C4 is drawn before P and C5 in one seeded stream, so removing C4
  sets there would have changed the BH-adjusted p of the planted sets and which P and C5 sets are drawn.)
- **Could it be steered?** I have not looked at any C4 rate.

### 4.3 Metabolomics C4 at size 20 (found while writing the first draft)

- **Declared:** metabolomics C4 at sizes 5, 10, 20; 5 sets per size. Reported; decides nothing.
- **L1's mechanics:** size 20 formed in 14 of 15 replicates (4 of 5 in cell B), identically under both models. Sizes 5
  and 10 formed in 15 of 15.
- **Arithmetic:** the model keeps every metabolomics feature (250 is fewer than the 1,000 it keeps per data set), so
  the pool is every clean complete feature in no block: L1's count of clean complete features per replicate (60 to
  83, `l1/clean-features.tsv`) less the 45 block members, that is 15 to 38 (mean 27, standard deviation about 5.6; 15
  in the replicate where the size did not form, 21 in the next smallest). The auditor's re-implementation: mean 27.6,
  standard deviation 5.6, at least 20 features in 0.93 of replicates. So the size is not one that cannot form: it sits
  just above the 90 % line and was flagged under the pre-L1 reading of the rule, per cell at five replicates.
- **Options:** (a) the reported rate over sizes 5 and 10. (b) Keep size 20 under the per-class reading.
- **Proposed: (a)**, in the scorer.
- **What it changes:** Nothing the simulation draws or tests changes: the sets are drawn as the scripts committed
  before L1 draw them (`simulate-k184.R` keeps C4 at its declared sizes), so every p-value, BH-adjusted p, power
  figure and K-E3 count is what the pre-L1 code gives. The scorer leaves the size out of the reported C4 rate and
  reports, for each size left out, how many replicates it formed in.
- **Could it be steered?** I have not looked at any C4 rate.

### 4.4 The feature-level rate on the main factor in metabolomics

- **Declared:** the feature-level rate on F_main "over the clean complete features the model did not select".
  Reported; decides nothing.
- **L1's mechanics:** 0 features in 15 of 15 replicates, under every arm (`l1/clean-features.tsv`).
- **Arithmetic:** the model keeps every metabolomics feature (250 against 1,000), so there is no feature it did not
  select. Empty by construction in every replicate of every cell.
- **Options:** (a) state it as not defined for metabolomics and report it for transcriptomics and proteomics only. (b)
  For metabolomics take all clean complete features, selected or not: these are features the factor was learned from,
  so it is a different quantity, and it would be a measure added after L1.
- **Proposed: (a).** The scorer marks those rows `defined = FALSE`; the simulation is unchanged.
- **What it changes:** one reported cell of the results is stated as not defined. The feature-level rates on F_ind
  and F_ret in metabolomics are unaffected. No criterion.

## 5. The scripts before a freeze, and the lifecycle extension L1b

**On the package branch, not run before L1b** (no R on the laptop; a scorer cannot be tested first without a run, so
L1b is what exercises this code before S; after the commit that carries this version the two scripts are not touched
again before S, or L1b is repeated):

1. `simulate-k184.R`: the error table has one row per set size (`nominal` added to its aggregation); every replicate
   writes its four package versions (K184-R01; before, only replicate 1 did, and S starts at 6). Nothing drawn or
   tested depends on either.
2. `score-k184.R`, the KIT-184 part: the sizes summed before every estimand; the scorer-side redefinitions of section
   4 as data (`left_out`, `c5_proteomics`, `feature_rate_undefined`); every check carries the composition it is read
   in, and one read in two compositions holds when both hold; `error-rates-by-size.tsv`; the determinism comparison
   (twelve tables, read as text, the counted jobs' rows of the same cell and replicate, both ordered by every column,
   `identical()`, `determinism.tsv`, one summary line); a summary line for R6; and three switches for L1b only
   (`K184_BLIND`: the summary is not printed; `K184_DETERMINISM_SELFTEST`: the comparison again with one value of one
   table altered, which must read FALSE; `K184_L1_REFERENCE`: cell A's rates with every size summed and nothing left
   out, against `l1/error-rates.tsv`: TRUE or FALSE, and, counts only, whether the numbers of sets and of replicates
   agree in every row and in how many rows T differs).

**On the lane branch** (`k184-evidence.yml`, `k184.R`): the blind L1b as one job, S's job list of section 3, the
`score` job's timeout and 90-day retention. `K184_BLIND` alone would not make L1b blind, because the scorer still
writes its rate tables into its output folder: the blindness is the lane's, which uploads one staging folder filled by
name, runs no `score` job and writes nothing to the step summary.

**L1b (lifecycle extension, blind; inside L1's ceiling, 52.3 job-minutes left; not a VAL-54 re-run).** One job on the
evidence lane, timeout 50 minutes, the scripts at the commit to be frozen and the lane's `k184.R` their byte copy. (1)
Cell A replicates 1 to 5, then cell D replicates 1 to 5, dealt as in S, with `K184_DEADLINE` at the start of the
simulate step plus 210 seconds. (2) Then, in the same job and without a deadline, cell A replicates 1 and 2 again into
a folder whose path contains `k184-det`. (3) Then `score-k184.R` in the same job with `K184_SCORE=ruv`,
`K184_BLIND=1`, `K184_DETERMINISM_SELFTEST=1` and `K184_L1_REFERENCE` set to the branch's `l1/error-rates.tsv`. (4)
Then the scale test: the worker tables copied with replicate numbers shifted until they hold about 2,000 replicates,
scored with `K184_BLIND=1` under `timeout 25m`, its output folder deleted; kept are its wall seconds, peak resident
memory and exit status. (5) Uploaded is only a staging folder filled by name: `plan.txt`, `runner.txt`,
`script-sha256.txt`, `status.tsv`, `completion.tsv`, `not-completed.tsv`, `formed.tsv`, `pools.tsv`, `versions.tsv`,
`determinism.tsv`, `determinism-selftest.tsv`, `l1-reference-check.tsv`, the scale test's line, the workers' logs, and
`tables.tsv` (for every worker and scorer table: name, rows, column names, SHA-256); the step prints the staged file
names. No `score` job runs and nothing is written to the step summary. Passes when: every cell A replicate completed;
at least one cell D replicate is `not_started` and counted; `ke1-checks` and `ke5-checks` have 34 rows and
`ke4-checks` 4 rows per deciding cell (cell B present with missing values); `error-rates-by-size` exists;
`versions.tsv` has four rows; `determinism.tsv` 12 of 12 TRUE and the altered copy FALSE; `l1-reference-check` TRUE;
the scale test ends within 25 minutes and below 12 GiB. A mechanical failure is fixed and L1b sent again under VAL-54
(at most three, inside L1's ceiling).

L1b is judged from its staged files alone, against that list. Nothing else of it is opened, downloaded or printed.

## 6. Register texts against the code, and two additions to L1b's output list

Settled by re-check 1 (its section (c)):

1. L1b's `not_started` replicate: the first text asked for two cell D replicates, which the lane would start together.
   Now five cell D replicates, dealt as in S (1, 2, 3, 4 to 5), with the deadline 210 seconds after the simulate step
   starts; the acceptance is at least one cell D replicate not started and counted, and every cell A replicate
   completed.
2. "Eleven tables" is twelve (R04's wording, used above).
3. "Each worker runs A or B before D" is R05's sentence about the lane's processes (used above).
4. The first L1b text could not be blind (the determinism replicates in a second job, the scorer's folder and the
   workers' error tables uploaded whole, a `score` job writing the summary): R02's one-job text replaces it.

Two additions of mine to R02's fixed list, both counts or failure text only, stated here before L1b is sent:

- `tables.tsv` has one more column, the rows per cell (for example `A=34;B=34`), for the tables whose first column is
  the cell. R02's pass list asks for row counts "per deciding cell", which a total alone does not show.
- If the scorer exits with an error (in its first run or in the scale test), the job prints the lines of the scorer's
  standard error from `Error` on (at most twelve): R's error message and call. Without it a failure of code that has
  never run could not be diagnosed inside the ceiling. The scorer prints no result on standard error; in a run that
  succeeds nothing of it is printed.

One change of mine to the scorer beyond R02's refinement, needed for R02's pass list: the check tables now list every
check with its composition whether or not the data filled it, so a deciding cell without data (cell B in L1b) shows 34
rows with missing values instead of 33.

Not contradicted, checked against the code: the register's seed formulas, the 33 checks, the order in which classes
are drawn (C0, C1, C2, C3, C4, P, then C5), one BH family per call.

Left open: the auditor's formation figures (metabolomics C5 at 93 %, the candidate-share window) are not recomputed by
me; whether S's runners match L1's processor is recorded per job and not controlled.
