# KIT-197 outcome

**Outcome O2 (power).** On 100 new draws per setting the higher false share among the proteins only the count rule
keeps is a power and composition effect: pooled E_comp 0.0131 (one-sided 97.5 % lower limit 0.0125, present), while
unchanged proteins among them reach p < 0.05 at the rate of those kept under both rules (ρ_O 1.017, 95 % upper limit
1.034, within the 1.10 margin: absent). The part is not specific to value dependence (present under random
missingness too, at about half the size). The pooled reading has one exception: unchanged added proteins holding a
value at or below the filter's intensity threshold (2.1 % of them pooled) reach p < 0.05 2.5 times as often (see
Inferred). Details and every estimate: `results.md`.

Ticket `tickets/kit/KIT-197.md` revision 1.2 (records), 1.3 (the implementation audit's text repairs I01 to I06,
`AUDIT-KIT-197-impl-20261124.md`; no new run, no table changed); branch `kit-197-filter-added-proteins` from
`v0.9-staging` at 19997323. Package code unchanged.

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

- The "biased group mean" named in VAL-111 shrinks the fold change of changed proteins in the added set toward 0
  (attenuation −0.021 against −0.006 log2) and does not move unchanged proteins' fold changes (m6_NA ≈ 0). It is
  small. The power gap (0.61 against 0.80) is also present with random missingness (0.63 against 0.78), where there is
  no shrinkage (attenuation 0.0003), so most of it comes with the missing value itself and with which proteins have
  one (not separated). The added sets differ between the arms (5.8).
- Values at or below the filter's threshold shift the fold change as derived (m6_ℓ 0.12) and carry a calibration
  excess: unchanged added proteins holding one reach p < 0.05 2.5 times as often as those kept under both rules [2.4,
  2.7]. They are 2.1 % of the added unchanged proteins in the pool, and 15 % (B vd5), 26 % (C vd5) and 36 % (D vd5) at
  the lowest missing share. Channels B2 and B3 are not separated (the design does not separate them).
- B vd5 reads "both" (ρ_O 1.17, lower limit 1.05), carried by those proteins (channel share 1.30 [0.88, 2.55]). It is
  seen also on S2's draws (ρ_O 1.10 [1.00, 1.20]). It is named in the sentence; heterogeneity is reported, not tested
  (5.8).

## Acceptance criteria

| # | Criterion | Met | Evidence |
|---|---|---|---|
| A1 | no package code changes | yes | `git diff --stat 19997323..branch`: only `docs/qa/workflow-execution/kit/kit-197/`, `NEWS.md` (one clause in KIT-187's bullet), `docs/methods-decision-ledger.md` (one row) |
| A2 | L1 after D0 to D5, on the scripts as declared, complete, blind; E0 to E6 pass; G2 applied | yes | `lifecycle.md`, `l1/checks.tsv`, run 45316713443; G2: 57 min planned, no split |
| A3 | declaration alone after L1, before S_a; hashes match | yes | commit 4786112c (one file); provenance of both jobs |
| A4 | S_a complete, G1 passes; S_b complete; determinism identical | yes | `runs.tsv`, `sa/results/g1.tsv`, `sb/rows/determinism.tsv`, both `checks.tsv` (2,400 of 2,400) |
| A5 | every estimand per setting with R, estimate, SEs, limits; descriptive decomposition; 12 labels; pooled label; attribution; outcome | yes | `sb/results/` (`settings.tsv`, `pooled.tsv`, `labels.tsv`, `r5-mcar.tsv`, `strata.tsv`, `m5.tsv`, `outcome.tsv`), `sa/results/descriptive.tsv` |
| A6 | outcome.md and results.md | yes, after the implementation audit's repairs | these files: observed and inferred separated, the threshold-value exception stated (I01), the power attribution corrected (I03), "For the change report" (I06); `results.md`: the Methods sentence as replaced (I02), M5 as worded and signed; the generated template stays in `sb/results/sentence.md` |
| A7 | NEWS clause and ledger row | yes, after the implementation audit's repairs | `NEWS.md` KIT-187 bullet, re-wrapped with the exception clause (I01 (c), I05); ledger row `VAL-149` after VAL-121 with I01 (d) and I04 (a) to (d) |
| A8 | core-hours per run | yes | `runs.tsv`: 10.54 of 25.33 for the counted runs; with the post-hoc R1-01 run (0.75) 11.29 |

## Departures (mechanics only)

- E6 kills the interrupted run with processx's `kill_tree()` once one task is complete, instead of a fixed `timeout`.
- L1 also exercises the determinism pass and the descriptive mode's G1 stop, which otherwise run first in S_a and S_b.
- The branch is `kit-197-filter-added-proteins` (the coordinator's name), not ticket 3's `kit-197-added-fdp-cause`.
- The template's "more | not more" unchanged proteins is decided by the 95 % interval of π0_Δ − π0_O (declared); it
  reads "more" (0.804 against 0.802) while the null-share part is −0.0008 (`results.md`).

## Post-hoc aggregate by ℓ (implementation audit R1-01; the owner's ruling of 2026-11-24)

- Purpose: put the intervals of the ℓ = 1 and ℓ = 0 aggregates (over the declared R7 strata, weights = tested
  unchanged proteins per stratum) in a committed table, computed by this ticket's own script `ell-aggregate-k197.R`
  from the committed rows, with the declared resampling (the scoring script's own definitions, evaluated unchanged).
- Post hoc: written after the results were seen. Reported, not deciding: no label, outcome or threshold depends on it,
  and the script asserts no expected value except its self-check (the committed pooled ρ_O estimate, limits and
  jackknife SE re-derived to a relative 1e-9, else it stops).
- Run: one request on the runner (`k197-ell-20e17c20`), binding 20e17c20 (the script as added in a935e8e7 plus one
  line creating its output directory, made before any run); ceiling 30 job-minutes (each script 15 minutes); outputs
  `audit/ell-aggregate-sb.tsv`, `audit/ell-aggregate-sa.tsv`. A failure in the mechanics follows VAL-54; a self-check
  failure is a result and is reported, not re-run.
- Observed: run 45342463572 (AMD EPYC 7763, reference BLAS), scripts job 11.2 minutes (0.75 core-hours; each script 3
  s), both self-checks passed at the first attempt. ℓ = 1: 2.535 [2.365, 2.709] (S_a 2.470 [2.288, 2.651]); ℓ = 0:
  0.984 [0.964, 1.004] (S_a 0.981 [0.960, 1.002]); random missingness ℓ = 1 1.553 [1.395, 1.713], ℓ = 0 1.012 [0.996,
  1.028]; ℓ = 1 share 0.0212 (S_a 0.0206); [2.5 %, 97.5 %] bootstrap limits. These equal the implementation audit's
  figures to the digits it printed.

## Not checked

- Real data; mean-variance trends, other missingness mechanisms, correlated proteins (5.8).
- The separate sizes of channels B2 (variance selection) and B3 (unequal group sizes).
- Whether the per-setting heterogeneity (B vd5, A vd5) is more than chance: not tested by design.
- G1 on the 9V45 processor S2 used (it held on an EPYC 7763).
- `fitFDistUnequalDF1` being the branch taken in every fit (not instrumented).
- The aggregate over strata of ρ_O for ℓ = 1 and ℓ = 0 is in no table the counted jobs wrote; `ell-aggregate-k197.R`
  (this ticket's script, run on the runner, 45342463572) computes it with its intervals from the committed rows
  (`audit/ell-aggregate-sb.tsv`, `audit/ell-aggregate-sa.tsv`), after re-deriving the committed pooled limits. It is a
  post-hoc aggregate over the declared strata, computed after the results were seen: reported, not deciding.

## With the owner (not decided here)

- Whether the package changes how the filter and the test treat values at or below the threshold, and whether the
  differential-analysis table gets a column marking the proteins kept with such a value (audit OD3).
- How the paper words the claim, for example the pooled result with the exception and the 4 x 3, 5 % setting named
  (audit OD1, OD4).

## For the change report

- What (observed in simulation): KIT-187 (VAL-71, VAL-108) changed the missing-value filter's count rule so that, with
  groups of 3, a group passes with one value missing. This ticket changes no code. It measures why the proteins kept
  only by the new rule have a higher false share among their calls. In simulation, pooled, about 6.4 % of them are
  false against 4.9 %, while the whole list stays at about 5 %. The cause is mainly lower power among them. The
  exception is the few whose compared groups hold a value at or below the filter's intensity threshold, which the
  filter counts as missing but the test uses: when unchanged they reach p < 0.05 2.5 times as often.
- Who is affected (inferred from the code and the simulation, not measured on users' data): proteomics analyses with
  groups of 3 or 4 and missing values, from KIT-187 on.
- Code: the threshold and the count of values at or below it as missing are at
  `R/func_prot_qc_filtering_helpers.R:229-233` and `:243-247` (19997323, unchanged at 3d4a585d). The test uses the
  masked matrix's finite values as they are (`R/func_prot_da_model_stats.R:306-313`).
- Results changed: none (records only). What users are told changes: the NEWS clause in KIT-187's bullet.
- Evidence: runs 45316713443 (lifecycle), 45319257326 (S2's draws), 45327303208 (new draws); `sb/results/pooled.tsv`,
  `labels.tsv`, `strata.tsv`. Settings: KIT-185's generator with no unwanted variation, cells 2 x 3, 4 x 3, 5 x 3 and
  3 x 4, 5, 10 and 25 % value-dependent and random missingness, 100 new draws per setting.
- Limits: simulated data only, from one generator, with equal group sizes; real data, unequal groups and mean-variance
  trends are not covered.


## Merge record (coordinator)

Merged into `v0.9-staging` on 2026-11-24 from 07b4cd2f after the implementation audit and its re-checks 1 and 2 (register `AUDIT-KIT-197-impl-20261124`): I01 to I06 and R1-01 closed. The ledger row is VAL-149. The owner's decisions on this result are VAL-142: the intervals on the threshold-value figures come from this ticket's own script run on a runner (run 45342463572); the threshold-value inconsistency is a recorded candidate, not fixed here; the paper's claim names the exception; the NEWS entry and the Methods sentence carry the exception, final wording at the owner's wording pass. No package code changed with this merge; NEWS gained one clause.
