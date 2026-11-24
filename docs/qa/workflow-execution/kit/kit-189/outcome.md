# KIT-189 outcome

Ticket revision 1.3; declaration `declaration.md` (commit 8deb1649, alone). Branch `kit-189-qvalue-vs-bh` from `v0.9-staging` at 3ff60190448150fa13e06d624378cad1481d3820. No package code changes.

**Outcome O1: the shipped q-value fails in 9 of its 36 deciding cells, Benjamini-Hochberg over the kept terms (BH) fails in none, and BH over the whole size-filtered family (BHf) controls in all 36. By rule F the recommended default is BHf. The q-value stays offered with its measured failure stated (VAL-86).**

## What was run

| Step | Runs | Commit | Result |
|---|---|---|---|
| L1 blind lifecycle | 45275836889; E7 instrument run 45281701515 | 426cebb0; 9f78d703 | E0 to E4, E6, E7 pass; G1 passed; G2 re-split (`lifecycle.md`), budget raised by the owner (VAL-101) |
| Declaration | — | 8deb1649 | committed alone; script hash equal in every S and D provenance file |
| S | 24 jobs (`runs.tsv`) | 8deb1649 | 22,800 replicates (45,600 route calls), every one written; 0 crashed, 0 empty calls; no timeout, no continuation, no instrument rerun |
| D (E5) | 45291725277 | 8deb1649 | the 280 rows of replicates 1 to 10 of cells 2 and 11 identical to S's without `part` and `seconds` |
| SC | 45311934798 | 8039b4e4 (rows committed) | once; outputs in `results/` |

Core-hours: 183.8 in all (`runs.tsv`), against 238.7 expected and 446.6 ceiling (VAL-101); S used 178.6.

## Observed

### Integrity (`results/checks-score.tsv`, `checks.tsv`)

- 319,200 rows; each replicate read once and complete; one list per replicate under both propagation values.
- E2's per-call checks hold on all 45,600 calls (margins, BH, the constant q / p.adjust, the chain reproducing every q-value, the phyper p-value, BHf calls within BH calls, BH calls within q-value calls). No list has every q-value NA. E5 passes.

### The 108 classifications (`results/classifications.tsv`; line 0.075)

Any-call rate on random lists: estimate (one-sided 95 % upper limit), classification.

| Cell | Propagation | Q | BH | BHf |
|---|---|---|---|---|
| B1 E. coli full, L 25 | none | 0.179 (0.193) **fails** | 0.072 (0.082) not determined | 0.020 (0.025) controls |
| | parent terms | 0.230 (0.246) **fails** | 0.051 (0.059) controls | 0.016 (0.021) controls |
| B1, L 126 | none | 0.160 (0.174) **fails** | 0.038 (0.045) controls | 0.016 (0.021) controls |
| | parent terms | 0.058 (0.067) controls | 0.026 (0.033) controls | 0.017 (0.023) controls |
| B1, L 400 | none | 0.038 (0.045) controls | 0.028 (0.035) controls | 0.025 (0.031) controls |
| | parent terms | 0.020 (0.026) controls | 0.017 (0.023) controls | 0.015 (0.020) controls |
| B2 E. coli 2,000 subset, L 25 | none | 0.085 (0.096) not determined | 0.056 (0.065) controls | 0.017 (0.022) controls |
| | parent terms | 0.167 (0.181) **fails** | 0.032 (0.039) controls | 0.017 (0.023) controls |
| B2, L 126 | none | 0.058 (0.067) controls | 0.025 (0.031) controls | 0.022 (0.028) controls |
| | parent terms | 0.031 (0.038) controls | 0.021 (0.027) controls | 0.018 (0.024) controls |
| B2, L 400 | none | 0.017 (0.022) controls | 0.016 (0.021) controls | 0.016 (0.021) controls |
| | parent terms | 0.021 (0.027) controls | 0.021 (0.027) controls | 0.021 (0.027) controls |
| B3 yeast, L 25 | none | 0.101 (0.113) **fails** | 0.083 (0.093) not determined | 0.010 (0.015) controls |
| | parent terms | 0.048 (0.057) controls | 0.048 (0.056) controls | 0.011 (0.015) controls |
| B3, L 126 | none | 0.429 (0.448) **fails** | 0.042 (0.050) controls | 0.012 (0.017) controls |
| | parent terms | 0.078 (0.088) not determined | 0.029 (0.036) controls | 0.015 (0.020) controls |
| B3, L 400 | none | 0.055 (0.064) controls | 0.024 (0.030) controls | 0.019 (0.025) controls |
| | parent terms | 0.025 (0.032) controls | 0.018 (0.024) controls | 0.014 (0.019) controls |

Mean false discovery proportion on planted lists of 126: estimate (one-sided upper limit), classification.

| Cell | Propagation | Q | BH | BHf |
|---|---|---|---|---|
| B1 small | none | 0.120 (0.142) **fails** | 0.026 (0.038) controls | 0.017 (0.027) controls |
| | parent terms | 0.050 (0.059) controls | 0.019 (0.025) controls | 0.013 (0.018) controls |
| B1 mid | none; parent | 0.053; 0.021 controls | 0.013; 0.008 controls | 0.006; 0.004 controls |
| B1 large | none; parent | 0.012; 0.003 controls | 0.002; 0.002 controls | 0.001; 0.001 controls |
| B2 small | none; parent | 0.033; 0.026 controls | 0.014; 0.019 controls | 0.010; 0.014 controls |
| B2 mid | none; parent | 0.024; 0.008 controls | 0.009; 0.004 controls | 0.005; 0.003 controls |
| B2 large | none; parent | 0.005; 0.003 controls | 0.003; 0.002 controls | 0.002; 0.001 controls |
| B3 small | none | 0.357 (0.388) **fails** | 0.019 (0.026) controls | 0.005 (0.008) controls |
| | parent terms | 0.061 (0.069) controls | 0.017 (0.021) controls | 0.007 (0.009) controls |
| B3 mid | none | 0.221 (0.241) **fails** | 0.014 (0.017) controls | 0.005 (0.008) controls |
| | parent terms | 0.038 (0.043) controls | 0.013 (0.016) controls | 0.006 (0.007) controls |
| B3 large | none; parent | 0.066 (0.074); 0.008 controls | 0.003; 0.002 controls | 0.001; 0.001 controls |

Totals: Q 9 fail, 25 control, 2 not determined; BH 0 fail, 34 control, 2 not determined (B1 and B3, lists of 25, no parent terms; upper limits 0.082 and 0.093, the largest rates not excluded); BHf 36 control.

Under the shipped propagation (`"is_a_part_of_regulates"`) the q-value fails in 2 cells (lists of 25 on the full E. coli proteome and on the subset; Bonferroni lower limits 0.199 and 0.140). Its other 7 failures are without parent terms (`"none"`, the value earlier records replay with, VAL-77).

### Power and the large-event rule (`results/paired.tsv`)

- Planted term called by Q minus by BH, paired: at most 0.039 in a deciding planted cell; the highest one-sided lower limit at 0.05 / 18 is 0.0005. **No material power gain.**
- Lists on which Q calls 10 or more terms and BH fewer than 10: up to 0.237 of random lists (yeast, 126, no parent terms; upper limit 0.253), and 0.151 with parent terms (E. coli full, 25; upper limit 0.165). The large-event condition (upper limit below 0.01 in every null cell) fails in 6 of 18 cells. On one random yeast list of 126 with parent terms, Q called 1,144 terms more than BH.
- On no list did Q call fewer terms than BH: 2.1 item 2 holds on every list.

### Mechanism (`results/mechanism.tsv`)

- pi0_hat (q / p.adjust) collapses: its median per null cell is 0.14 to 1.00; lists with pi0_hat below 0.5 occur in 15 of 18 null cells (up to 1,684 of 2,000: yeast, 126, no parent terms; smallest value 0.001).
- The within-list counterfactual (A-Qf: the same chain over the whole size-filtered family, count-0 terms entered at p = 1): its median pi0_hat is 1.00 in every null cell, and on every list where the kept-family pi0_hat is below 0.5 the whole-family one is at least 0.9, except 3 lists in one cell (E. coli subset, 400, no parent terms). Dropping the count-0 terms is what collapses pi0, as 2.1 inferred.
- Stage of enrichit's chain: at L = 25 the second stage (bootstrap at lambda 0.05; no kept p-value at or above 0.95) gives 14 to 100 % of the q-values per cell; at L = 126 and 400 the first stage (the smoother) gives most. A-Q's any-call rate is far higher on stage-1 lists (for example 0.50 against 0.04 on yeast, 126, no parent terms).

## Inferred

- The q-value's failure is a property of the route (the kept family; read in 2.1 and confirmed by the counterfactual), not of a user's data. It is largest without parent terms and for small lists, and remains with parent terms at L = 25.
- Restricting BH to the kept terms loosens control at L = 25 (BH's two "not determined" cells, estimates 0.072 and 0.083, against BHf's 0.020 and 0.010). BHf costs some power: planted term called in 74.8 % of tested lists at the lowest, against 78.5 % for BH and 80.5 % for Q.

## Recommendation put to the owner (O1, rule F; VAL-86)

- **Default: Benjamini-Hochberg over the whole size-filtered family (BHf).** A package change: the family is formed by the package, not by enrichit.
- **Offered: the q-value**, with its measured failure stated.
- BH over the kept terms (the existing `p.adjust` column) is not named by O1; its text is given so the owner may also offer it (VAL-86 asks for "the q-value and Benjamini-Hochberg").
- The report and help say which flag a run used; the application's notice ("adjusts p-values with BH") is corrected with the change.

Each option's trade-off text T: `results/options.tsv`, reproduced in `results.md`.

## Acceptance criteria

| # | State | Evidence |
|---|---|---|
| A1 | met | `git diff --stat 3ff60190..kit-189-qvalue-vs-bh` lists only `docs/qa/workflow-execution/kit/kit-189/` and `docs/methods-decision-ledger.md` |
| A2 | met | `lifecycle.md`, `checks.tsv`; runs 45275836889, 45281701515 |
| A3 | met | commit 8deb1649 holds `declaration.md` alone; every provenance file's `script_sha256` equals the declared `sim-k189.R` hash; SC ran `score-k189.R` unchanged since the declaration |
| A4 | met | 22,800 of 22,800 replicates; no crash; E5 pass (`runs.tsv`, `rows/`, `rows-d/`) |
| A5 | met | `results/`; outcome O1 and its recommendation |
| A6 | met | this file; `results.md` |
| A7 | met | ledger row `VAL-<merge>` |
| A8 | met | `runs.tsv`: 183.8 core-hours |

## Not checked

- Real lists, down lists, cutoffs other than 0.05, other organisms or annotation sources, the two-relation propagation set.
- The g:Profiler route; the integration's legacy enricher fallback (its `p.adjust` is BH over kept pathways, so BH's two "not determined" cells bear on it; not measured).
- Whether DOSE's historical code also adjusted over terms with hits only.
- The literature cited in 2.1 (Tarone 1990, Gilbert 2005, Benjamini and Yekutieli 2001) was not read; `results.md` does not cite it.
- `results.md`'s Methods wording is proposed; the owner edits it (VAL-41).
