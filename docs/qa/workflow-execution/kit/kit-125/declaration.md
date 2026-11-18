# KIT-125: declaration before any run (copy of the ticket's revision 1.4, 2026-11-18)

The ticket file is not in the repository; this is its revision 1.4 as written before the lifecycle test and the confirming runs. Revisions 1.1 to 1.3 (the design, the definition of a contradiction, the user's decisions on the survey) are quoted in the ticket and summarised in `outcome.md`.

## Declaration (fixed before any run)

**Purpose.** Do the verdicts of KIT-124's simulations about three designs hold on real noise? One null split study per omic.

**What is run (revision 1.1, M-01).** The arm code of KIT-124's simulators, on the tree the simulations used: work tree `kit/wt-k124-sim` at `0b6fcb6d` (read-only), `R_LIBS = kit/scratch-k124-lib : core/browser-smoke.y5ordS/library`, `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, R 4.6.1, as `kit/scripts-k124/resume-20261118.sh` sets them. The package kernels of the recommended design are not merged (`kit-126` classic proteomics, `kit-127` lipid and metab; DIA_LIMPA and transcriptomics are being written on `kit-129-k-dia` and `kit-129-trx`), so for every omic the recommended design is the simulator's own arm-7 code path, copied into the KIT-125 script; today's design and no RUV are the package's kernels as the simulators call them. This makes the real-data figures comparable with the simulated ones; it does not test the unmerged kernels.

| Design | Simulator arm | Classic proteomics (S1) | DIA_LIMPA (S2) | Transcriptomics (S3) | Metabolomics (S4) |
|---|---|---|---|---|---|
| No RUV | 1 | `runTestsContrasts` | the route's `dpcDE` test | `trxDaFit`, `trxDaApply` | the default path's limma kernel |
| Today's design | 3 | `findBestNegCtrlPercentage`, `ruvIII_C_Varying`, test | the sweep, `computeRuvCovariates`, test on the corrected EList (correction mode) | `trxRuvGrid`, `trxRuvs`, `trxDaContrast` | `runPerAssayRuvOptimization`, `ruvIII_C_Varying`, `runMetabolitesDA` |
| Today's default mode (DIA_LIMPA only) | 4 | | the same sweep, W in the model (covariates mode, the route's shipped default) | | |
| Recommended | 7 | all proteins as controls, groups as replicate sets, k = E, W in the model | the same | `RUVSeq::RUVs`, all genes, groups as replicate sets, k = E, W in the model | all features as controls, groups as replicate sets, k = E, W in the model |

E is `sva::num.sv(method = "be", B = 20)` as each simulator calls it; E = 0 makes the recommended design no RUV for that split. Arm 4 is added for DIA_LIMPA because arm 3 there is the route's correction mode and the route's shipped default is the covariates mode; both are "today".

**Data sets** (accessions opened at their source on 2026-11-18; the manifest with URL, licence, sha256, size and date is `data-manifest.tsv`).

| Omic | Accession | What | Units of one condition | Batches | Licence and what it allows |
|---|---|---|---|---|---|
| Classic proteomics | Zenodo 10.5281/zenodo.17607826 (source data PXD023175) | FragPipe `combined_protein.tsv`, MaxLFQ, plasma, timsTOF | 30 healthy volunteers (biological) | none (one plate) | CC BY 4.0: use and redistribution with credit |
| DIA_LIMPA | PRIDE PXD024043 (Brunner et al. 2022) | DIA-NN 1.8 main report, single HeLa cells | 203 untreated cells (label UB; biological) | two culture batches, measured on two dates | CC0 |
| Transcriptomics | ENA PRJEB5348; counts from github.com/bartongroup/profDGE48 at commit `5ff79149` (Schurch et al. 2016) | gene counts, yeast wild type | 48 biological replicates, 42 after the authors' exclusion list | none recorded | repository MIT; ENA data open |
| Metabolomics | MetaboLights MTBLS404 (Thevenot et al. 2015, "sacurine") | 120 identified urine metabolites, LC-HRMS negative mode | 184 adult volunteers (biological) | ne1 (90), ne2 (94) | EMBL-EBI terms of use (no explicit open licence): used for computing and reporting; the table is not copied into the repository |

The stub asks for one data set per omic, "lipidomics or metabolomics": metabolomics is taken (MTBLS404: biological units, two recorded batches, no licence restriction on reuse). Lipidomics ST003456 (CC BY-NC-ND, no batch recorded), MaxQuant PXD002854 and TMT Zenodo 14767905 (technical replicates only), and DIA with IQ on PXD024043, all accepted in revision 1.3, are not run in this pass; they are named in the outcome as not checked. Nothing derived from any data set goes into the repository beyond per-split aggregates (counts, rates) and the public sample names of each split.

**Preparation** (fixed here, blind to the split results; it stands where a simulator generates its data).

- Classic proteomics: the `MaxLFQ Intensity` columns of the samples named `HC_*`; 0 is missing; log2; no further normalisation. Pool: the HC samples with at least half the median number of quantified proteins among the 30 (a rule against failed injections; the counts run from 73 to 1,084). Per split, S1's filter: a protein is kept when at least half of each false group's samples are observed.
- DIA_LIMPA: the report's rows with `Q.Value`, `PG.Q.Value` and, where the columns exist, `Global.Q.Value` and `Global.PG.Q.Value` at most 0.01 (the shipped thresholds); `Precursor.Normalised` above 0, log2; protein = `Protein.Group`, precursor = `Precursor.Id`. Pool: the runs labelled UB with at least 500 precursors (the shipped per-sample cutoff). Per split: precursors observed in at least 2 runs of each false group (the shipped `min_reps_per_group = 2`, `min_groups = 2`, as counts), then the route's DPC quantification on the split's 12 runs, engine `fastlimpadpc` (the shipped engine, as S2), no protein normalisation (shipped: none).
- Transcriptomics: the 48 count files joined; the counting tool's five summary rows dropped; pool = 42 replicates. Per split the package's own filter and TMM, as S3.
- Metabolomics: the MAF's columns of material type "biological material"; 0 or empty is missing; log2; one assay; no normalisation. Per split, S4's presence rule (observed at least once in each false group).

**Split unit and design (revision 1.1, M-02).** The unit is the biological unit (a volunteer, a cell, a culture replicate); no data set here has technical replicates. Per split, 12 units are drawn from the pool and put into two false groups of 6, the group size of the simulations' centre cell. Two split types:

- random: 12 units drawn from the pool;
- batch-aligned, where two batches are recorded (DIA_LIMPA, metabolomics): group A = 4 units of batch 1 and 2 of batch 2, group B = 2 and 4. The correlation of batch with the false group is 1/3 by layout, the simulations' unbalanced-batch layout and close to their cell C5 (0.3). A full alignment (A = batch 1, B = batch 2) is not run: the batch difference would then be the group difference and no design could separate them.

**Number of splits and seeds.** 100 distinct splits per omic and split type (the stub's number; within revision 1.1's ceiling of 200; at revision 1.2's threshold of 100 for "informative"), numbered 1 to 100. Seed of a split = 5,000,000 + 100,000 x omic index (classic 1, DIA_LIMPA 2, transcriptomics 3, metabolomics 4) + 10,000 x type index (random 0, batch 1) + split number; the units are drawn under that seed; arm a runs under seed x 100 + a; E uses its simulator's stream. Two equal partitions in one run stop the run before any arm. The lifecycle test uses splits 9001 to 9003.

**Measures (KIT-124's rules 1 and 7, revisions 1.1 and 1.2).** Per split and design: the share of tested features with p under 0.05 (null rate) and under 0.01, and whether the omic's rule calls anything (q under 0.05). Per omic, split type and design: the mean null rate, its range over the splits, and the share of splits with any call. The splits share samples, so no confidence bound is put over them (revision 1.1); the range and the 5th and 95th percentiles of the per-split null rate are the "interval" shown. A refused split carries its no-RUV result; a crashed split has no measures and is counted.

**Limits (revision 1.2).** Mean null rate at most 0.075. Share of splits with any call: at most twice no-RUV's share on the same splits, or at most 0.10. Stated limit of the method: a null split measures false calls, not power; one data set is one unit.

**Which designs "passed" KIT-124's rules 1 and 7** (from `kit-124/scored/s{1,2,3,4}-core/k124-cells.tsv`, full subset; "failed" = a fail of rule 1 in a judged core cell or of rule 7; fixed here before any run):

| Omic | No RUV (arm 1) | Today (arm 3) | Recommended (arm 7) |
|---|---|---|---|
| Classic proteomics | passed | failed (rule 1 in six cells; rule 7: 0.206) | passed (rule 7 inconclusive: 0.061, upper bound 0.078) |
| DIA_LIMPA | failed (rule 1 in C5 only, 0.116; rule 7 passed, 0.047) | failed (arm 3 and arm 4: rule 1 in nearly every cell; rule 7: 0.576, 0.564) | passed (rule 7: 0.044) |
| Transcriptomics | passed | failed (rule 1 in C5 only, 0.179; rule 7 inconclusive: 0.078) | passed (rule 7 inconclusive: 0.074) |
| Lipid, metab | passed | failed (rule 1 in every judged cell) | passed (rule 7 has no valid check: C4 is void in S4 by its scoring unit) |

**Decision each outcome leads to (revision 1.2, applied by the letter).**

- A passed design above a limit on real data: **contradiction**. Stop: KIT-129 does not proceed; to the user.
- A failed design within both limits on every split type of its data set: **contradiction**. The same stop.
- Anything else: **consistent**. KIT-129 proceeds as far as this ticket is concerned.
- Fewer than 100 scored splits for an omic and split type (crashes, or a run that stopped): **not informative**; reported as not checked.

Foreseen before the run and written here so that it is not an explanation after the fact: two "failed" designs failed rule 1 only in cell C5 (the factor tied to the groups) and were within the limits in the simulations' own global-null cell C4: transcriptomics' today's design (C4 null rate 0.056) and DIA_LIMPA's no RUV (C4 0.052). A random split has no factor tied to the groups in expectation, so the simulations themselves predict these two within the limits on random splits. If that happens, the letter of revision 1.2 says "contradiction" for transcriptomics (no batch-aligned split exists there) and the outcome reports it as such, prominently, with the matching simulation cell beside it; whether it stops KIT-129 is the user's call. The simulation cell matched to each split type (C4 for random; C5 and C4 for batch-aligned) is printed beside every real figure.

**Budget and attempts.**

- Lifecycle test first (AGENTS.md): 3 splits per omic and split type, splits 9001 to 9003, through completion, with one deliberate failure (`K125_FAIL_AT=9002` stops the recommended design's arm in that split) to prove that a failure is recorded as a crash with its message and that the run goes on and finishes. Laptop, one R process, `nice -n 19`. Debugging of the scripts happens here, on these splits only; split numbers 1 to 100 are not run before the confirming run.
- Sizing from the lifecycle's seconds per split. An omic whose 100 splits per type are sized at one core-hour or less runs on the laptop: one R process, `nice -n 19`, a `timeout` of 5,400 s per omic (one and a half times the size limit, for a loaded machine). An omic sized above one core-hour is not run on the laptop; it is prepared for GitHub-hosted runners (the public files are downloaded by the runner from the source; nothing is committed) and reported.
- Attempts: one per omic. A run that fails for an instrument reason (a script error, the timeout, a crash share that leaves fewer than 100 scored splits) is recorded and that omic stops; no retry, no change of preparation, seeds or limits after the run.
- No memory figure is taken.


## Bounded check of the unwanted variation's strength (added 2026-11-18, after the DIA_LIMPA result and before any of the figures below was computed)

**Why.** DIA_LIMPA's no RUV was scored "failed" in simulation (cell C5, null rate 0.116) and is within both limits on PXD024043, also on the batch-aligned splits: a contradiction by the rule above, recorded and accepted by the user as not bearing on the default (decision M21, `kit-124/design-review-20261118.md` section 6). M21 asks for one bounded check of one reading: **the labelled unwanted variation of this data set (the two culture batches) is weaker than the unwanted factor S2 generates.** This section fixes how that is measured. It changes no rule, verdict or figure above.

**Matrices.** Twelve samples each, quantified as the null splits and the simulations were (the route's DPC quantification, engine `fastlimpadpc`, on the pinned tree `0b6fcb6d`, the libraries above).

- Real: the batch-aligned splits 1 to 10 of the run above (their units redrawn under their seeds and checked against `runs/dia_limpa/splits.tsv`), with `k125_dia_limpa.R`'s own preparation, filter and `quantify`. Each holds 6 cells of each culture batch. The factor is the batch indicator.
- Simulated: S2's own `simulateData` and `quantify`, read from `kit-124/sim/s2_limpa.R` at the pinned tree and not edited, cells C0 and C5, runs 1 to 10 under the simulator's seeds (the first ten runs KIT-124 scored). The factor is the true score vector the generator returns. Only the 800 proteins without a group effect are measured, since the real false groups have none.

**Measures, per matrix.**

1. Per protein, the share of the variance of its 12 quantified values explained by the factor (R-squared of the regression on the factor): median, quartiles, and the shares of proteins above 0.10 and above 0.30.
2. Pooled: the factor's share of the total protein-centred sum of squares of the matrix.
3. The floor of 1 and 2 with 12 samples: the same two measures with the factor's values permuted over the samples once (seed 100 x the matrix's seed + 98).
4. Principal components of the protein-centred matrix: the share of variance of each of the first five, and each one's R-squared with the factor.
5. Confounding: the correlation of the factor with the B-against-A indicator (real: 1/3 by layout; C5: the realised value; C0: 0 in expectation).

Reported side by side as the median over the ten matrices with the range.

**Reading, fixed here.** Decided on measure 1's median and on measure 2:

- **supported** (the real batch is weaker): on both measures every real matrix is below every simulated run of C5, and of C0;
- **refuted** (the real batch is not weaker): on either measure the median over the real matrices is at or above the median over C5's runs. This is stated prominently;
- **undetermined**: anything else (the ranges overlap).

Measures 3 to 5 are shown for the reader and decide nothing.

**Limits, stated before the result.** One data set. Batch is the only unwanted factor of the real data that is labelled; unlabelled unwanted variation is not measured, so "weaker" can only ever mean "the labelled batch is weaker". The real factor is a balanced two-level indicator and the simulated one a continuous score. Twelve samples per matrix, so a factor with no effect still explains about 0.09 of a protein's variance on average (measure 3 gives the floor). The check measures strength; it does not show that strength is what moves no RUV's null rate.

**Budget.** Laptop, one R process, `nice -n 19`, `timeout` 1,800 s, one attempt (`k125_strength.R`, output under `runs/strength/`). A script error before the table is written is an instrument failure: one correction of the script, recorded, with no change to this section.
