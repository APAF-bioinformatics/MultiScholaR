# KIT-124 revision 1.1: independent statistical review, 2026-11-16 (condensed register)

Reviewer: a separate agent that did not write the ticket; read-only; it ran no simulation. **Verdict: not ready to run.** As written the likeliest output for every omic is "inconclusive, so RUV off by default" whatever the truth (R-01, R-04, R-06), and a pass would not show calibration in the sense RUV is meant to deliver (R-02, R-03). Every repair is an edit to the ticket. Not yet applied: revision 1.2 follows W0.

## Blocks the decision

| ID | Finding | Fix |
|---|---|---|
| R-01 | A good design cannot pass rule 2 at the planned precision: at SE 0.0075 and limit 0.02 a design with zero true excess passes one cell with probability 0.76 and 16 cells with 0.012. | SE at most 0.0045 for rule 2, with R-02. |
| R-02 | Rule 2 compares against no-RUV's realised false discovery proportion: the noisiest reference (SD 0.14 against 0.025 for the oracle design) and the wrong one. A calibrated design can fail where no-RUV is conservative; an inflated one can pass where no-RUV is itself inflated (factor tied to group, batch left out). | Bound the mean false discovery proportion absolutely (upper bound at most 0.075); report paired differences against the true-factor arm, and the batch arm in the batch cell. |
| R-03 | The rules test calibration averaged over experiments, and "benefit" asks for no benefit. With the factor independent of group, no-RUV is calibrated on average while 5 to 20% of its runs exceed 10% false discoveries; a design that does nothing passes. | A per-run rule (share of runs above 10%, against the true-factor arm plus a margin); rule 3 becomes a gain over no-RUV. Efron 2007; Leek and Storey 2008. |
| R-04 | Rule 4's 0.02 sits on the arithmetic cost of two extra covariates at 6 per group (0.016 to 0.023 for k = 2; 0.005 to 0.010 for k = 1), and was set after the no-factor run. | The owner sets the tolerated cost knowingly, or "k = 0 must be reachable" is the requirement. |
| R-05 | Arm 13 names no estimator; k is undefined for most arms. With all features as controls the package's rule returns 1. | One named estimator that can return 0 (parallel analysis on group-residualised data); a k-source column. |
| R-06 | The fallback (no RUV) is not held to the rules and one hard cell vetoes everything; no check that a cell is passable; inconclusive undefined below 1,000 runs. | Judge arm 1 by the same rules; core cells (must pass) and stress cells (state conditions); void cells the true-factor arm does not pass. |

## Fix before running

| ID | Finding | Fix |
|---|---|---|
| R-07 | Rules 5 and 6 need exact binomial bounds and a number of runs (0 failures in 200 has an upper bound of 0.018). Refusals and crashes counted alike. | Clopper-Pearson; at least 1,000 data sets for the null split; refusals counted apart and scored as no-RUV results. |
| R-08 | Repetition rule ambiguous (pilot at the centre only; SDs range 0.01 to 0.23 across cells); the compute ceiling contradicts it. | Pilot every cell or fix 1,000; cut arms before runs. |
| R-09 | The null split is 3 against 3 at the centre, has no no-RUV reference and does not test control selection with true effects, a factor tied to group, or power. | State the limits; add a 0% true-effects cell at full group size; compare with arm 1. |
| R-10 | The simulators are RUV's own model and under-specified (S2, S3 one line each; factor strength; how the correlation with group is made; symmetric effects independent of loadings favour RUV; "50 known nulls" is an oracle, not spike-ins). | A parameter table per simulator; cells with asymmetric effects and with loadings tied to abundance; arm 9 relabelled "oracle controls" plus a degraded supplied set. |
| R-11 | In S2 the factor must be added at precursor level before missingness is drawn, or limpa's route is flattered. `simProteinDataSet` defaults to 25 proteins and ignores `sigma`. | Add to `other$E.complete`, redraw missingness with the returned DPC coefficients. |
| R-12 | One setting at a time misses the interactions that matter (3 per group with a correlated factor or 3 factors; no factor with 50% missing). The audit's F-10 asked for corners; the remediation recorded it as applied without them. | 4 to 6 corner cells as stress cells. |
| R-13 | `RUV4`, `RUVinv` and `sva` take no missing values; their p-value and k source are unnamed. | Score all arms on the complete subset as well; name `p.ebayes`. |
| R-14 | The reference-usage column mislabels two things: for RUVs the vignette passes all genes; the RUV-III and RUV-III-C papers tested on corrected values with technical replicates. The tie-break "closest to reference usage" is undefined. | W0 confirms at source; the preference order per omic is written before any run. |
| R-15 | "Higher power" at the nominal cut rewards the more liberal design. | Compare power at matched achieved false discovery proportion. |
| R-16 | The test's settings differ by omic and are not fixed in the simulators. By reading: the proteomics default passes a TREAT cutoff of 0, takes the TREAT branch and so appears to drop trend and robust. | Settings per simulator as the default path produces them, pinned by a test. The proteomics observation is recorded, not acted on. |

## Minor

R-17 Bradley's criterion used one-sided and at 0.05 only. R-18 "simulated and real data" against real data out of scope: narrow the claim or set parameters from aggregates. R-19 true effects are the first rows while control selection breaks ties by row order; seed scheme; one random stream per arm. R-20 sources to add (Gagnon-Bartsch and Speed 2012; Leek and Storey 2007; Gerard and Stephens; De Livera; Kim 2021; Korthauer 2019; Siepe 2024).

## Found sound

The structure after Morris et al.; one-sided bounds in every cell (conservative against a false pass); the arm contrasts that isolate the test form; the comparators; all features as controls as a documented usage; every range in the ticket's exploratory table recomputed from the files; the code facts of the implementations table.

## Could not check

Bradley 1978 (paywalled); Morris, Weber, Poulos, Molania read through summaries; S2 to S4 do not exist yet.
