# KIT-211 outcome

Ticket `tickets/kit/KIT-211.md` r1.4 (measurement; no package code changed). Branch `kit-211-threshold-filter-test`, cut from
v0.9-staging 894b9cc821a70b03baab1e1fa47b3d47440c7ed5 (bound commit; package tree equal to 4d588ba3's). Every request's `package_commit`:
cfcde043d1daadb44e18b18f8e6ce83cf7914662 (the ten counted jobs and L1's successor state) and 656b08c65cd3cac19af0aab3725178b7b9ffb7fd (the scoring; it adds only the
returned rows); both have an empty `git diff --stat 894b9cc8 <commit> -- R/ src/ inst/ tests/ DESCRIPTION NAMESPACE`.

## What ran

- G0 (`l1/g0.txt`): the three protein filter files unchanged since 19997323; no stop. DA-file changes (KIT-200, KIT-201, KIT-220, KIT-249) scope E1; `R/func_lipid_qc.R` (KIT-202) scopes E9.
- Declaration `declaration.md` committed before L1 (bffa7859; hashes in 15829c8f), unchanged afterwards. No VAL-54 section was needed.
- L1: run 45713532325, all checks pass at the first attempt (`lifecycle.md`, `l1/`).
- Ten counted jobs (all complete, exit 0, every task done; `rows/*/range.tsv`; determinism D: 5 of 5 identical): 45723161360 (A 201-300, with D), 45723175176 (A 301-400), 45723180223 (B 201-300),
  45723175202 (B 301-400), 45723178712 (C 201-300), 45723204611 (C 301-400), 45723185718 (D 201-300), 45723192772 (D 301-400), 45723199616 (Av 201-400), 45723199520 (Ag 201-400).
- Scoring: run 45730910292, `score-k211.R` with no environment variables; 63 s. checks.tsv read first: complete 36,000 of 36,000; ranges incomplete 0; F1 0; F2 0; t_r monotone failures 0; BH fallbacks 0;
  determinism 5 of 5; bootstrap_B 10,000 (checks.tsv and outcome.tsv).
- No instrument defect, no re-run, no continuation, no S′ or SC′ reserve used.

## Observed (separated from inferred; the declared labels only)

- Outcome row **K2**. Positive control passes (Φ_S 1.4897, one-sided 97.5 % lower limit 1.4301 > 1.20).
- Classes: A **C** (primary "removes"; power cost "not determined": Δpow 0.0049408, one-sided limits [0.0047333, 0.0051567] straddle 0.005; guard "not raised"); B **E**; B′ **E**. Limits are in `results.md` (the one-sided 99.17 % φ limits: A [0.7636, 0.9183], B [0.9872, 1.0490], B′ [1.1042, 1.1791]).
- Falsifiers F3, F4, F5 (B and B′) did not fire; F1 and F2 did not fire (`results/falsifiers.tsv`, `checks.tsv`).
- Decision table, per-setting labels, near-threshold figures, whole-list figures, the reported cells Av and Ag: `results/decision-table.md`, `results.md`, `results/settings.tsv`, `results/pooled.tsv`.
- Per-setting labels differ from the pooled label in 22 deciding setting-candidate labels; they are named in `results.md`.
- **Per-setting qualification (RA-01).** The pooled values are the prespecified analysis, and they are within the 1.20 margin (B 1.02, B′ 1.14). In every one of the ten deciding settings, however, Φ_B and Φ_B′ are above 1 (B 1.04 to 1.38; B′ 1.11 to 1.48), with an excess remaining (one-sided 0.05/72 labels) in 3 settings for B and 7 for B′. The pooled values sit below the settings' because the arms' subgroups are distributed differently across the settings (a Simpson-type pooling effect of ratios of sums): B's and B′'s subgroups are concentrated in the 5 % and 10 % settings of the 4 x 3 and 5 x 3 cells (B vd5 and C vd5 hold 41 % of B's subgroup against 23 % of the no-filter subgroup), where every arm's rate is lowest, while the no-filter subgroup is spread almost evenly, the vd25 settings holding 30 % of it against 13 % of B's. A is above 1 in the 2 x 3 cell (1.35 to 1.58), where the shipped filter's own excess is not resolved (subgroups of 0.4 to 5 proteins per data set). Never read the pooled label as "B leaves no excess".
- **A's costs.** A's whole-list power change is 0.0049 (one-sided limits [0.0047, 0.0052] around the 0.005 margin, not determined) and its attenuation near the threshold is -0.144 against -0.045 shipped (about three times as large).

## Definitions the declaration left implicit (DR-03)

- ρ^{ℓ0}_X is computed on the whole kept set H_X against ℓ0 = 0 (R-c's form), not on Δ_X against O.
- The attenuation in R-h's band is taken over the band's tested planted proteins; its power is over every planted protein of the band (filter losses count against power).
- The jackknife SE is KIT-197's unstratified formula (delete one replicate of one cell with all its masks and arms), reported beside the bootstrap, which governs the labels.

## Acceptance criteria

| # | Criterion | Evidence |
|---|---|---|
| A1 | No package code changes | `git diff --stat 894b9cc8..<branch>` lists only `kit-211/` and `docs/methods-decision-ledger.md` (checked at the tip) |
| A2 | G0 before L1; L1 once, blind, complete; E0 to E10 pass; G2 applied | `l1/`, `lifecycle.md`, run 45713532325 |
| A3 | Declaration committed before any counted job; nothing declared changed | commits bffa7859, 15829c8f; `declaration.md` |
| A4 | Ten counted jobs complete; F2 and E10 in every replicate-mask; D identical | `runs.tsv`, `rows/*/range.tsv`, `checks.tsv` (F1 0, F2 0, determinism 5 of 5) |
| A5 | Every estimand per arm with R, bootstrap and jackknife SEs and limits; PC; each candidate's states and class; per-setting labels; one outcome | `results/` |
| A6 | Decision table and this file | `results.md`, `results/decision-table.md` |
| A7 | Ledger row | `docs/methods-decision-ledger.md` (VAL-<merge>) |
| A8 | Core-hours per run against section 11 | `runs.tsv`: counted jobs 426.7 job-minutes = 28.45 core-hours; L1 23.4 min = 1.56; scoring about 15 min = about 1.0 (script 63 s; not queried); about 465 job-minutes, about 31.0 core-hours in all. The counted-job minutes are request-created to run-last-updated times, an upper bound of the jobs. Against 38.51 restated expected and 53.28 declared; ceiling 103.49 / 105.60 |

At outcome time, `git diff --stat 894b9cc8 7bbc344e` over G0's files (staging head 7bbc344e, as the results audit ran it): only `R/func_lipid_qc.R` (KIT-217; the lipid helper, used here only by the lifecycle check E9) and `R/func_prot_da_model_methods.R` (KIT-204; "Not tested" in count tables) changed. Both are off the measured path, so every figure measured at 894b9cc8 also describes staging's code at 7bbc344e.

## With the owner

The outcome table's K2 asks the owner which candidate (B or B′; A is C): the choice between B and B′ is whether the percentile control keeps a meaning (6.9), not a measured difference; the tie-break is not applied between B and B′. The flag column (SD8), the Methods sentence, and the proteomics "groups required" default are decided with the direction, with Ag's reported figures in `results.md`. No direction is chosen here.

## For the change report

The proteomics protein missing-value filter counts a finite value at or below its intensity threshold (the 1st percentile of every finite value of the table) as missing, while the limma test uses that value (R/func_prot_qc_filtering_helpers.R:224-233, :242-247; R/func_prot_da_model_methods.R:287, :403-412; R/func_prot_da_model_stats.R at 91e3cdcb; lines moved with KIT-200, KIT-201). Proteomics DIA without limpa, LFQ and TMT meet it at every setting, and more proteins are exposed since KIT-187. KIT-197's 2.5-fold figure mixes this with an effect of how the subgroup is defined; this measurement splits it: the shipped ratio 2.60 = 1.49 (the filter's factor, Φ_S, against the unfiltered analysis; one-sided 97.5 % lower limit 1.43) x 1.75 (the ratio without a filter); on complete data, proteins holding such a value reach p < 0.05 1.22 times as often as proteins holding no such value (the part no filter rule can remove). Pooled over ten settings; in every individual setting Φ_B and Φ_B′ stay above 1 (B 1.04 to 1.38, B′ 1.11 to 1.48; excess remaining in 3 and 7 settings), the pooled values being lower because the subgroups are distributed differently across settings. This ticket changes nothing in the package. What each direction would change, with figures, settings and limits: `results.md` (A: whole-list power 0.0049 lower, near-threshold attenuation larger; B and B′: whole-list power 0.0037 and 0.0032 higher, about 24 and 20 more proteins kept per data set; Φ 1.02 and 1.14). Limits of the claim: one generator whose variance does not depend on intensity (Av reported), equal groups, TREAT with trend and robust settings, no normalisation or RUV between filter and test, the default percentile and minimums, strict mode not run.

## Not checked

- Real data; the intensity-dependent variance beyond the one reported cell Av; unequal groups; other percentiles; TREAT at a non-zero cutoff; the declared designs; S's and A's filter at 1 group.
- Whether any later QC step (normalisation, RUV-III-C, duplicate removal) changes a threshold-level value before the DA (the implementation's first gate).
- Bootstrap limits are not hand-checked (R's RNG is not emulated); only degenerate and exact quantities were (E4).
- The simulation's `K211_PARTIAL=1` path and a continuation on a fresh runner were not exercised (no range was incomplete).
- Job minutes in `runs.tsv` are request-created to run-last-updated (the API rate limit prevented per-job queries), so they are upper bounds of the jobs.
