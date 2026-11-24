# KIT-184 simulation: declaration, DRAFT

**This is a draft. It is not the declaration, it fixes nothing, and nothing may be counted against it.** The counted
run S has not started and does not start from this file. Written 2026-11-24 on the coordinator's ruling after the
lifecycle run L1 stopped at a stop rule (`lifecycle.md`): the draft holds the declaration as I would freeze it
(section 2), the changes L1 forces with their arithmetic (section 3), what I could not resolve (section 4) and what
would have to change in the scripts before a freeze (section 5). An independent audit reads it first. The declaration
itself is then `declaration.md`, committed alone, before S.

Ticket `tickets/kit/KIT-184.md` revision 1.5, sections 6.2 to 6.5 and 11; ledger VAL-53, VAL-58, VAL-62, VAL-94 (KIT-183's
outcome O1), VAL-97 (3).

## 1. What I had seen of L1 when I wrote this

The changes of section 3 alter what the simulation measures, and I propose them after seeing L1's results. So that a
reader can judge whether a proposal could have been steered, this is everything of L1 (run 45333818847, 5 replicates of
each of cells A, B and D) that I looked at, from the session's own record of my commands.

**Rejection rates and other results I saw** (all from 5 replicates per cell; all in `l1/` on the branch):

- The scorer's summary (`l1/results-summary.md`): how many of K-E1's 33 checks passed in cells A and B (17 and 18) with
  the largest upper limit; K-E3's counts; K-E4's rate and upper limit per data set in cells A and B (in cell B all three
  data sets fail it: 0.080, 0.107, 0.200); K-E5's pass counts (19 and 21 of 33); the outcome rows.
- `l1/error-rates.tsv`: class C5 only, test `camera`, arms 1, 2, 3 and 3o, every cell, on F_main, F_ind and F_ret,
  **pooled over set sizes** (the table in `lifecycle.md`, and the same for F_ret). For proteomics C5 that means I saw,
  for arm 3: 0.000, 0.107, 0.020 on F_main and 0.093, 0.120, 0.147 on F_ind (cells A, B, D).
- `l1/feature-rates.tsv`: cell A on F_ind and F_ret, every arm; arm 3 on F_main in every cell. Four rows of one
  worker's `feature.tsv` with their counts of significant features.
- `l1/power.tsv`: arm 3, `camera` and `camera_floor`, sets of 10, F_main, cells A and B (what [PF] and [PC] read).
- `l1/numbers.tsv`: the nine text numbers as 5 replicates would give them, with where each was found.
- `l1/caught.tsv`, `l1/identity.tsv`, `l1/exactness.tsv`, `l1/floor.tsv`, `l1/completion.tsv`, seconds per replicate.

**From this I wrote, in `lifecycle.md` and in my report, an indication that K-E4 may fail in cell B** (arm 3 calls C5
sets on the main factor in 8 to 20 %, the oracle arm in 0 to 1 %).

**What I did not look at:** any rate split by set size, for any class; any rate of C4, C4b, C1-out or C3; the rows of
`ke1-checks.tsv` and `ke5-checks.tsv` (only the summary's pass counts); `paired-differences.tsv`;
`operating-characteristic.tsv`; any per-replicate rate.

**Mechanics I used for section 3** (counts only): `l1/pools.tsv` and the workers' `formed.tsv` (sets wanted and formed
per replicate, class and size); the `features` column of the workers' `feature.tsv` (how many clean complete features
entered the feature-level rate per replicate); the generator's constants in `simulate-k184.R`. The script did not write
pool sizes: where a pool size is given below it is arithmetic from the generator's constants or from those counts, and
says so.

**Lifecycle replicates are excluded from S.** L1 ran replicates 1 to 5 of each cell, and the generator's seed depends
only on the cell and the replicate number. S therefore runs replicates 6 to 805 of cells A and B and 6 to 405 of cell D
(800, 800 and 400 replicates, none of them seen), and the determinism job replicates 6 to 15 of cell A. No L1 output
enters any count, rate or text number. (The ticket words the determinism job as "replicates 1 to 10": see section 4.)

## 2. The declaration as I would freeze it

Sizes marked † are changed by section 3; everything else is the ticket's, as the scripts already hold it.

### Vehicle, cells and jobs

- MOFA2 evidence lane `.github/workflows/k184-evidence.yml` on branch `run-mofa-evidence-k184` (the lane's
  `.github/mofa-smoke/k184.R` is a byte copy of `simulate-k184.R`, checked by the job): R 4.6.1, the KIT-134 snapshot
  library, MultiScholaR at the package commit named in the lane's `k184-params.json`, MOFA2 1.22.0, ruv 0.9.7.1,
  RUVIIIC 1.0.19, limma 3.68.4 (as L1 found them, `l1/versions.tsv`), one thread per worker
  (`OPENBLAS_NUM_THREADS = OMP_NUM_THREADS = MKL_NUM_THREADS = 1`), four workers per job.
- Scripts: `simulate-k184.R` and `score-k184.R` at that commit, their SHA-256 written here at the freeze (they change
  once more before it, section 5). KIT-183's two files are not edited.

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
  cell A replicates 6 to 15 again, timeout 25 minutes. Why no job holds cell D alone: under four busy workers a cell D
  replicate took 120 to 133 s in L1; 32 of them in one worker would start the last one after about 66 of the 68.4
  minutes available. Mixed, a job's last replicate starts after about 60 minutes.
- A replicate in which MOFA2, the correction or a test fails is `did_not_complete`. A cell with more than 5 % of its
  declared replicates not completed (failed or not started) is "not decided".

### Generator (per replicate r of a cell)

One stream, `set.seed(20261840 + offset + r)` with offset 0 (A), 100,000 (B), 300,000 (D). Everything of KIT-183's
declaration (`kit-183/declaration.md`, "Generator") is unchanged: features (4,000; 1,800; 250 on 204 identifiers),
abundances, noise, Z1 on 10 % of loading units and Z2 on 5 % of the others (none in metabolomics), value-dependent
missing values in proteomics and metabolomics, the metabolomics identifiers carrying two or three features. Added:

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
within the call. `gsea` seed 20261121; `permutation` 1,000 reorderings, seed 20261185 + r. `camera_floor` and
`permutation` are reported only (VAL-97 (3)). Under every arm that runs both, p_floor at least p_camera is checked set
by set.

### Classes (members are identifiers of the test's universe; drawn per replicate, model, data set and factor)

| Class | Members | Sets per replicate | Sizes: transcriptomics; proteomics; metabolomics |
|---|---|---|---|
| C0 | uniform from the universe | 10 per size | 10, 30, 100; 10, 30, 100; 5, 10, 20 |
| C1 | each formed null block, by rho_block | 2; 1; 1 per rho_block and size | 10, 30, 100; 10, 30, 100; 5, 10 |
| C2 | the lowest-abundance fifth (metabolomics: third) of the clean complete features in no block | 5 per size | 10, 30; 10, 30; 5 |
| C5 | complete features that are not candidates, zero loading on Z1 and Z2, positive loading on U1 | 5 per size | 10, 30, 100; **10, 30 †**; 5 |
| P-up, P-down | half with a positive (negative) Z1 loading, half clean complete in no block | 5 per size and direction | 10, 30, 100; 10, 30; 5 |
| C4 | clean complete features in no block that the model selected | 5 per size | corrected model **10, 30 †**, default model **none †**; 10, 30; **5, 10 †** |
| C1-out, C4b | block members the model did not select (at least 5; not metabolomics); those it selected (at least 5) | per block | as the blocks |
| C3 | metabolomics identifiers carrying two or three zero-loading features | 5 per size | 5, 10 |

A set of a size forms in a replicate when its pool holds at least that many features; sets of one class and size are
drawn independently from the same pool.

### Estimands (as ticket 6.3)

- Error rate: s_r = the share of the replicate's sets of the class (all sizes together) with unadjusted p at most
  0.05; T = the mean of s_r over the replicates with a set; interval T plus or minus t(R - 1) sd(s_r) / sqrt(R).
- Paired differences T1 minus T3 and T2 minus T3 per class, data set and factor, on F_ind, F_ret and F_main.
- Feature-level rate per arm (1, 2, 3, 3o): the share of clean complete features with moderated-t p at most 0.05 for
  the factor, on F_ind, F_ret and F_main; on F_main over the clean complete features the model did not select
  (**not defined for metabolomics †**).
- Power: the share of P-up and P-down sets with BH-adjusted p at most 0.05 in the expected direction on F_main, per
  arm, test, size, direction and data set; direction accuracy with the scorer's own orientation.
- Exactness per data set: the status of the written factors, the largest s_(k+1)/s_1, the largest relative residual.
  The identity check per replicate, data set and factor (arm 3's `correlation` against the same model on the corrected
  values; below 1e-8). What the removed factors caught (smallest canonical correlation with the planted factors).
- Determinism: the determinism job's result tables identical to the same replicates in S.

### Criteria, outcomes and the numbers in the texts (as ticket 6.4 and 6.5; nothing added, no threshold changed)

- **K-E1** (decides in A and B): arm 3, `camera`, upper 95 % limit of T at most 0.075 on F_main and F_ind for C0, C1
  (each rho_block) and C2, and on F_ind for C5, in every data set: 33 checks per cell.
- **K-E2** (`camera_pr`, C0 and C2): decides nowhere under KIT-183's O1; run and reported.
- **K-E3**: direction accuracy of arm 3 at least 0.99 in cells A and B, decided when at least 50 planted sets were called.
- **K-E4** (declared, does not hold the merge): arm 3, `camera`, the same limit on C5 on F_main, every data set, A and B.
- **K-E5**: arm 2, `camera`, the same limit on K-E1's 33 checks.
- Outcomes R1 to R8 exactly as ticket 6.5. An unfavourable valid result is reported, not re-run.
- [R3], [I3], [R5], [I5], [R6], [I6], [R5m], [PF], [PC] as ticket 6.4 and 5.5 define them, as the scorer's `numbers.tsv`
  gives them.
- Reported, deciding nothing: arms 1, 2, 3o, D0, D1 on every class; C4, C4b, C1-out, C3; the paired differences; the
  feature-level rates; power; cell D; `camera_floor` and `permutation`; the exactness measures; the operating
  characteristic from the measured standard errors.

### Attempts

One counted run S. One correction S' only for a failure of S's instrumentation, limited to re-running the affected
jobs within 305 job-minutes (20.33 core-hours); if a whole repeat were needed, stop, to the owner (ticket 11, revision
1.5). An unfavourable valid result is recorded and its outcome applied.

### Dated section (after L1, before S), as it would read

- Lifecycle: run 45333818847 (15 of 15 replicates completed; parity with KIT-183's cell identical; 110 to 133 s per
  replicate under four workers). The stop rule of ticket 11 fired at 70 s; the coordinator accepted the re-split above
  and kept the ceiling (ticket 11, revision 1.5). The first push (45333717225) ran no job.
- Classes and rates redefined after L1's formation counts (section 3 of the draft, once audited): proteomics C5 sizes
  10, 30; C4 sizes for transcriptomics 10, 30 under the corrected model and none under the default model, for
  metabolomics 5, 10; the feature-level rate on F_main not defined for metabolomics. K-E1 keeps 33 checks per cell and
  its limit.
- S runs replicates 6 onward; L1's replicates enter nothing.

## 3. Changes forced by L1

The ticket's rule (6.2, last sentence; KIT-183's rule): "a class with fewer than 90 % of its declared sets formed is
redefined in a dated section of the declaration before the counted run." How it is applied here: the scorer committed
before L1 ran (`score-k184.R` at 21053957, `pool_share = 0.9`) flags a row per cell, model, data set, factor, class and
size, and I take a size as below the rule when it is flagged in any cell. With 5 replicates per cell that means "did
not form in at least one replicate of a cell". KIT-183's scorer also counts per cell, and its dated section drops a
size that formed 80 %. Section 4, Q1, gives the other reading and what it would change.

The coordinator named three points. Writing this I found a fourth (3.3), which `lifecycle.md` had missed and now states.

In every case below the remedy I propose is to remove what cannot form. I add no size, class, criterion, margin, cell
or arm, and change no threshold.

### 3.1 Proteomics C5 at size 100

- **Declared:** C5 in proteomics at sizes 10, 30 and 100, 5 sets per size (15 sets per replicate). The ticket's own
  indicative pool is "about 100" features (6.2, pools paragraph), written before any run.
- **L1's mechanics:** size 100 formed in 10 of 15 replicates (4, 3, 3 of 5 in cells A, B, D), identically under both
  models, since the pool does not depend on the model. Sizes 10 and 30 formed in 15 of 15.
- **Arithmetic** (generator constants; the pool size itself was not written): 1,800 features, 180 with a Z1 loading,
  90 with a Z2 loading, 1,530 with neither; half of those (765) are candidates; of the other 765 a share 0.4 x 0.5 = 0.2
  has a positive loading on U1, about 153; at the ticket's indicative completeness of 0.674 about 103 are complete
  (binomial standard deviation about 9.5). A pool of 103 plus or minus 9.5 holds 100 features in about 65 % of
  replicates; L1 found 67 %. Two more things follow from the same arithmetic: where size 100 does form, its five sets
  are five draws of 100 from about 100 to 115 features, so they are nearly one set counted five times; and as declared,
  s_r would be a share over 15 sets in about two thirds of replicates and over 10 in the rest.
- **Options:** (a) proteomics C5 at sizes 10 and 30 only. (b) Enlarge the pool by lowering the proteomics candidate
  share from 0.50: by the same arithmetic a share between about 0.43 and 0.45 would let both the blocks (420 complete
  candidates needed) and C5 at size 100 form in at least 90 % of replicates. This changes the generator for every
  proteomics class and arm (the share of proteomics features carrying unwanted variation rises from 36.8 % to about
  40 %), rests on a completeness figure L1 did not measure for these pools, and would need another lifecycle run to
  confirm. (c) Replace 100 by another size: a new size chosen after L1. (d) Keep the size and average over the
  replicates in which it forms: the rule forbids it, and the estimate would be conditional on the pool.
- **Proposed: (a).** It is the rule's own remedy and KIT-183's precedent, it chooses no new number, and it leaves the
  generator alone. Transcriptomics keeps C5 at size 100 (pool about 476 by the same arithmetic; formed 15 of 15).
- **What it changes:** the proteomics C5 rate becomes the share over 10 sets per replicate (five of 10, five of 30).
  That rate is read by one of K-E1's 33 checks (proteomics, F_ind: deciding), by K-E5's same check, by K-E4's
  proteomics check, and is a candidate for [R3], [R5], [R6] and [R5m]. The number of checks (33), the limit (0.075)
  and every other class are unchanged. With fewer sets per replicate the interval of that one check is, if anything,
  slightly wider.
- **Could it be steered?** I had seen proteomics C5's pooled rates (section 1), which in 5 replicates sit above the
  limit for arm 3 on F_ind and, in cell B, on F_main. I had not seen any rate by set size, so I do not know whether
  removing size 100 moves that rate up or down. The option that keeps size 100, (b), is the one that would change the
  data every arm sees.

### 3.2 Transcriptomics C4

- **Declared:** C4 ("as KIT-183 forms them", on the restated pool): clean complete features in no block that the model
  selected; transcriptomics sizes 10, 30, 100; 5 sets per size. Reported for every arm; decides nothing; no text number
  reads it.
- **L1's mechanics:** under the default model size 10 formed in 2 of 15 replicates and sizes 30 and 100 in none. Under
  the corrected model sizes 10 and 30 formed in 15 of 15 and size 100 in 14 of 15 (4 of 5 in cell B).
- **Arithmetic:** the model keeps 1,000 of the 4,000 transcriptomics features, those with the largest variance above
  the trend. About 1,877 features are clean (the 1,020 candidates and 0.36 of the other 2,380 zero-loading features; L1
  counted 1,841 to 1,926 per replicate), of which about 1,037 are in no block. Default model (values before RUV): about
  1,907 features carry at least one unwanted loading of 1 to 2 noise standard deviations, more than the 1,000 the model
  keeps, so clean features are almost never kept and the pool is nearly empty by construction. Corrected model: from
  L1's counts of clean complete features with and without those the model kept, the model kept 306 to 362 clean
  features in 14 replicates, of which about 0.55 (170 to 200) are in no block; in one replicate of cell B it kept 77
  (about 43 in no block), so size 100 did not form there.
- **Options:** (a) per model: corrected model sizes 10 and 30; default model no C4 for transcriptomics (arms D0 and D1
  then report no transcriptomics C4, and the declaration says why). (b) One definition for both models: every size is
  flagged somewhere, so no transcriptomics C4 at all. (c) Redefine C4's pool as KIT-183 had it before unwanted
  variation existed (zero loading on Z1 and Z2, whatever the unwanted loadings): under the default model that pool is
  mostly features that share unwanted variation, a different class under the same name. (d) Under the other reading of
  the rule (Q1) the corrected model would keep size 100.
- **Proposed: (a).** It removes only what did not form, and keeps what the corrected model can measure.
- **What it changes:** reported rows only: no transcriptomics C4 under arms D0 and D1; no size 100 in the
  transcriptomics C4 rate of arms 1, 2, 3, 3o. No criterion, no check, no text number.
- **Could it be steered?** I have not looked at any C4 rate.

### 3.3 Metabolomics C4 at size 20 (found while writing this draft)

- **Declared:** metabolomics C4 at sizes 5, 10, 20; 5 sets per size. Reported; decides nothing.
- **L1's mechanics:** size 20 formed in 14 of 15 replicates (4 of 5 in cell B), identically under both models. Sizes 5
  and 10 formed in 15 of 15.
- **Arithmetic:** the model keeps every metabolomics feature (250 is fewer than the 1,000 it keeps per data set), so
  the pool is every clean complete feature in no block: L1's count of clean complete features per replicate (60 to
  83) minus the 45 block members, that is 15 to 38 (mean 27, standard deviation about 5.6; 15 in the replicate where
  the size did not form, 21 in the next smallest). At that mean and spread about one replicate in eleven has fewer
  than 20.
- **Options:** (a) sizes 5 and 10. (b) Keep size 20 under the other reading of the rule (Q1): by the arithmetic it
  would then form in about 90 % of S's replicates, on the line.
- **Proposed: (a).**
- **What it changes:** reported rows only: no size 20 in the metabolomics C4 rate. No criterion, check or text number.
- **Could it be steered?** I have not looked at any C4 rate.

### 3.4 The feature-level rate on the main factor in metabolomics

- **Declared:** the feature-level rate on F_main "over the clean complete features the model did not select".
  Reported; decides nothing.
- **L1's mechanics:** 0 features in 15 of 15 replicates, under every arm.
- **Arithmetic:** the model keeps every metabolomics feature (250 against 1,000), so there is no feature it did not
  select. Empty by construction in every replicate of every cell.
- **Options:** (a) state it as not defined for metabolomics and report it for transcriptomics and proteomics only. (b)
  For metabolomics take all clean complete features, selected or not: these are features the factor was learned from,
  so it is a different quantity, and it would be a measure added after L1.
- **Proposed: (a).**
- **What it changes:** one reported cell of the results is stated as not defined. The feature-level rates on F_ind
  and F_ret in metabolomics are unaffected. No criterion.
- **Could it be steered?** There is no number to steer.

## 4. Not resolved, or needing a second pair of eyes

- **Q1: per cell, or over L1's 15 replicates?** I read the 90 % rule per cell, as the scorer was coded before L1 and
  as KIT-183 applied it. Read over all 15 replicates instead, two reported rows would stay as declared: transcriptomics
  C4 at size 100 under the corrected model and metabolomics C4 at size 20 (each 14 of 15, 93 %). Proteomics C5 at size
  100 (67 %) and transcriptomics C4 under the default model (0 to 13 %) fall under either reading. Nothing deciding
  depends on the answer.
- **Q2: the determinism job's replicate numbers.** The ticket says "replicates 1 to 10 of cell A". With L1's
  replicates excluded I propose 6 to 15. Engineering, for the coordinator to confirm.
- **Q3: transcriptomics C4 at size 100 under the corrected model failed for a different reason than the others.** Its
  pool is usually well above 100; it failed in one replicate in which the model kept few clean features. I infer, and
  have not shown, that this is a replicate in which the correction left unwanted variance in the corrected values. If
  so, a rate over the replicates where the size forms would lean towards replicates the correction served well, which
  is a further reason to drop the size. I did not look at that replicate's rates to check.
- **Not changed, though the arithmetic is tight:** nothing else. Every other class and size formed in 15 of 15.

## 5. What would change in the scripts before a freeze (not done)

1. `simulate-k184.R`: proteomics C5 sizes; C4 sizes per model and data set. Both scripts' SHA-256 are then re-pinned.
2. `score-k184.R`: **the determinism comparison is missing.** The KIT-184 part of the scorer leaves the determinism
   job's files out of the rates, as it should, but never compares them with S's (KIT-183's part does, for its own
   cells). Found while writing this draft; an instrument gap, to be closed before S.
3. The lane's job list (section 2) and `k184-params.json`.
4. Because 1 and 2 change code that L1 did not exercise, I propose a second lifecycle run before the freeze, on
   replicates L1 already saw so that nothing new is looked at: cell A replicates 1 to 4 and a determinism job on
   replicates 1 to 2, to show the changed classes form, the scorer writes its determinism table and a difference would
   be reported. About 27 job-minutes, inside what is left of L1's ceiling (75 less the 22.7 used). For the coordinator
   to decide.
