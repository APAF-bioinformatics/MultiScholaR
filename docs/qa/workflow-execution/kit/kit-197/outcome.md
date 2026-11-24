# KIT-197 outcome

**Outcome O2 (power).** On 100 new draws per setting the higher false share among the proteins only the count rule
keeps is a power and composition effect: pooled E_comp 0.0131 (one-sided 97.5 % lower limit 0.0125, present), while
unchanged proteins among them reach p < 0.05 at the rate of those kept under both rules (ρ_O 1.017, 95 % upper limit
1.034, within the 1.10 margin: absent). The part is not specific to value dependence (present under random
missingness too, at about half the size). Details and every estimate: `results.md`.

Ticket `tickets/kit/KIT-197.md` revision 1.2; branch `kit-197-filter-added-proteins` from `v0.9-staging` at
19997323. Package code unchanged.

## What was run

| Step | Run | Commit | What | Result |
|---|---|---|---|---|
| 0 | none | (ticket text) | the re-check's repairs K197-R1-01 to R1-06 applied with their exact texts; revision 1.2; closure table section 19 | done |
| scripts | none | 5a702b7b | `sim-k197.R`, `score-k197.R`, `l1-k197.R`, `e0-k197.py`, `e4hand.py` | E0 passed locally (one stated line differs) |
| L1 | 45316713443 | 5a702b7b | blind lifecycle E1 to E6, determinism and G1 mechanics | all pass at the first attempt (`lifecycle.md`) |
| declaration | none | 4786112c | `declaration.md` committed alone | hashes match what S_a and S_b ran |
| S_a | 45319257326 | 4786112c | r = 1 to 100 (S2's draws), descriptive scoring | G1 96 of 96; F1 0; descriptive tables |
| S_b | 45327303208 | 4786112c | r = 101 to 200 (new draws), determinism pass, deciding scoring | O2; F1 0; F4 0; determinism 5 of 5 |

No instrument defect, no rerun, no reserve used. Core-hours (job minutes × 4 / 60): L1 1.25, S_a 4.43, S_b 4.86;
total **10.54** against 10.67 expected and 25.33 ceiling (`runs.tsv`).

## Observed

- G1: every S2 figure (mean |O|, |H|, |Δ|, calls, power, fdp, `added_share_called`, `added_fdp`, and the replicate
  counts) reproduced to a relative 1e-9 in all 12 vd settings, on an AMD EPYC 7763 runner (S2 ran on a 9V45).
- Script hashes in both counted jobs' provenance: `sim-k197.R` 71001de9…c599, `score-k197.R` f5469479…e426, as
  declared; bound commit 4786112c in both.
- Pooled (ten deciding settings, vd): ρ_O 1.0169 [0.9968 one-sided 0.025; 1.0340 at 95 %]; E_comp 0.01313 [0.01248;
  0.01368]; label power; jackknife-normal label power. Mcar pool: ρ_O 1.0213 (present and absent: "non-zero, within
  the margin", so absent), E_comp 0.00646 present; label power. Attribution of P: not specific to value dependence.
- Per-setting labels (vd5 / vd10 / vd25): A not determined / power / power; B both / power / power; C power / power /
  power; D bias / neither / power (D vd5, vd10 not deciding). Mcar: A power ×3; B not determined / power / power; C
  not determined / not determined / power; D bias / bias / power. Bootstrap and jackknife labels agree in all 24.
- m6_NA pooled −0.0002 [−0.0020, 0.0015]; no setting's Bonferroni interval excludes 0 (F4 not triggered).
- S_a (descriptive, no labels): ρ_O 1.011 [0.991, 1.033], E_comp 0.0135 [0.0129, 0.0142], E_cal 0.0017
  [−0.0008, 0.0042]: the same picture on S2's own draws.
- t_r monotone in every replicate-mask (the counterfactual's premise); q-values by `qvalue` in all 4,800
  replicate-masks (no BH fallback).

## Inferred

- The "biased group mean" named in VAL-111 acts, in this generator, through power: changed proteins in the added set
  have their fold change shrunk toward 0 more (attenuation −0.021 against −0.006 log2) and are found less often
  (power 0.61 against 0.80); the missing value does not move unchanged proteins' fold changes (m6_NA ≈ 0).
- Values at or below the filter's threshold shift the fold change as derived (m6_ℓ 0.12) but do not show as a
  calibration excess in the pool; channels B2 and B3 are not separated (the design does not separate them).
- B vd5 reads "both" (ρ_O 1.17, lower limit 1.05): a null-side excess in one deciding setting, named in the sentence;
  heterogeneity is reported, not tested (5.8).

## Acceptance criteria

| # | Criterion | Met | Evidence |
|---|---|---|---|
| A1 | no package code changes | yes | `git diff --stat 19997323..branch`: only `docs/qa/workflow-execution/kit/kit-197/`, `NEWS.md` (one clause in KIT-187's bullet), `docs/methods-decision-ledger.md` (one row) |
| A2 | L1 after D0 to D5, on the scripts as declared, complete, blind; E0 to E6 pass; G2 applied | yes | `lifecycle.md`, `l1/checks.tsv`, run 45316713443; G2: 57 min planned, no split |
| A3 | declaration alone after L1, before S_a; hashes match | yes | commit 4786112c (one file); provenance of both jobs |
| A4 | S_a complete, G1 passes; S_b complete; determinism identical | yes | `runs.tsv`, `sa/results/g1.tsv`, `sb/rows/determinism.tsv`, both `checks.tsv` (2,400 of 2,400) |
| A5 | every estimand per setting with R, estimate, SEs, limits; descriptive decomposition; 12 labels; pooled label; attribution; outcome | yes | `sb/results/` (`settings.tsv`, `pooled.tsv`, `labels.tsv`, `r5-mcar.tsv`, `strata.tsv`, `m5.tsv`, `outcome.tsv`), `sa/results/descriptive.tsv` |
| A6 | outcome.md and results.md | yes | these files; M5 as worded and signed in `results.md` |
| A7 | NEWS clause and ledger row | yes | `NEWS.md` KIT-187 bullet; ledger row `VAL-<merge>` after VAL-121 |
| A8 | core-hours per run | yes | `runs.tsv`: 10.54 of 25.33 |

## Departures (mechanics only)

- E6 kills the interrupted run with processx's `kill_tree()` once one task is complete, instead of a fixed `timeout`.
- L1 also exercises the determinism pass and the descriptive mode's G1 stop, which otherwise run first in S_a and S_b.
- The branch is `kit-197-filter-added-proteins` (the coordinator's name), not ticket 3's `kit-197-added-fdp-cause`.
- The template's "more | not more" unchanged proteins is decided by the 95 % interval of π0_Δ − π0_O (declared); it
  reads "more" (0.804 against 0.802) while the null-share part is −0.0008 (`results.md`).

## Not checked

- Real data; mean-variance trends, other missingness mechanisms, correlated proteins (5.8).
- The separate sizes of channels B2 (variance selection) and B3 (unequal group sizes).
- Whether the per-setting heterogeneity (B vd5, A vd5) is more than chance: not tested by design.
- G1 on the 9V45 processor S2 used (it held on an EPYC 7763).
- `fitFDistUnequalDF1` being the branch taken in every fit (not instrumented).
