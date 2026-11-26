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
| A8 | Core-hours per run against section 11 | `runs.tsv`: ten counted jobs about 28.4 core-hours (request created to last update, an upper bound of the job), scoring about 1.0 (63 s script; request-to-update minutes in `runs.tsv`), L1 1.56 |

At outcome time `git diff <bound> <current staging head>` over G0's files: not run here (the staging head moved after the bound commit); every figure is "measured at 894b9cc8".

## With the owner

The outcome table's K2 asks the owner which candidate (B or B′; A is C): the choice between B and B′ is whether the percentile control keeps a meaning (6.9), not a measured difference; the tie-break is not applied between B and B′. The flag column (SD8), the Methods sentence, and the proteomics "groups required" default are decided with the direction, with Ag's reported figures in `results.md`. No direction is chosen here.

## For the change report

The proteomics protein missing-value filter counts a finite value at or below its intensity threshold (the 1st percentile of every finite value of the table) as missing, while the limma test uses that value (R/func_prot_qc_filtering_helpers.R:224-233, :242-247; R/func_prot_da_model_methods.R:287, :403-412; R/func_prot_da_model_stats.R at 91e3cdcb; lines moved with KIT-200, KIT-201). Proteomics DIA without limpa, LFQ and TMT meet it at every setting, and more proteins are exposed since KIT-187. KIT-197's 2.5-fold figure mixes this with an effect of how the subgroup is defined; this measurement splits it (Φ_S 1.49 against the unfiltered analysis; 1.75 with no filter; 1.22 on complete data; shipped 2.60 against proteins kept under both rules). This ticket changes nothing in the package. What each direction would change, with figures, settings and limits: `results.md` (A: whole-list power 0.0049 lower, near-threshold attenuation larger; B and B′: whole-list power 0.0037 and 0.0032 higher, about 24 and 20 more proteins kept per data set; Φ 1.02 and 1.14). Limits of the claim: one generator whose variance does not depend on intensity (Av reported), equal groups, TREAT with trend and robust settings, no normalisation or RUV between filter and test, the default percentile and minimums, strict mode not run.

## Not checked

- Real data; the intensity-dependent variance beyond the one reported cell Av; unequal groups; other percentiles; TREAT at a non-zero cutoff; the declared designs; S's and A's filter at 1 group.
- Whether any later QC step (normalisation, RUV-III-C, duplicate removal) changes a threshold-level value before the DA (the implementation's first gate).
- Bootstrap limits are not hand-checked (R's RNG is not emulated); only degenerate and exact quantities were (E4).
- The simulation's `K211_PARTIAL=1` path and a continuation on a fresh runner were not exercised (no range was incomplete).
- Job minutes in `runs.tsv` are request-created to run-last-updated (the API rate limit prevented per-job queries), so they are upper bounds of the jobs.
