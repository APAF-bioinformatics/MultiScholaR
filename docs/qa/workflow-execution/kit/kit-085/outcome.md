# KIT-085 outcome: the lipid and metabolite RUV sweeps impute once per assay

**The lipid and metabolite sweeps now impute the canonical-correlation input once per assay, as proteomics does, and the results are 1:1 with imputing at every percentage. This makes no corpus we have any faster, because none has missing values. Its value is on data with missing values, which no journey has exercised.**

User direction, 2026-11-05: "optimise it to the same standard as proteomics."

## The change

`prepareOmicsRuvCancorObject` (`R/func_metab_norm.R`) gives each sweep an object whose assays carry their NIPALS imputation, prepared once before the percentages run (f8327f97). Only the canonical-correlation plots read that object. The control selection and the correction keep the data with its missing values. No public signature changed, and `ruvCancor` itself is unchanged.

An assay the preparation cannot handle stays as it is, so `ruvCancor` handles and reports it as before.

## Evidence

- **New test:** `test-workflow-omics-ruv-cancor-once.R`, 16 pass on real sweeps of both omics. Before the change, 5 percentages over 2 assays imputed 10 times; now 2. Results are identical to the old path on 5% missing values.
- **Targeted RUV files:** 209 pass, 0 fail, including the new test.
- **Dependents:** 42 files, 7,199 pass, 0 fail, 0 error, 39 skip (`dependents.txt`).
- **`R CMD check`:** at baseline, 3 WARNINGs and 3 NOTEs (`check.log`).
- **Screen journeys on `kit_085` = RQRIvK, against `kit_082`'s oracles:**

| Member | Decision | Exact | Time | Whole slice | Application |
|---|---|---|---|---|---|
| `kit-journey-085l` | accepted | yes, 208 files | 355 s | 1.83 GiB | 1.02 GiB |
| `kit-journey-085m` | accepted | yes, 199 files | 339 s | 1.86 GiB | 0.94 GiB |

Memory equals `kit_082`'s. Both logs show no NIPALS run and 40 "No NAs", so these journeys show no regression and do not exercise the new path. The 1:1 test does.

The large journeys were not rerun: they would not reach the new path either.

## Where RUV's time goes

From the worker logs of the large journeys on `kit_084` (`ruv-timings.txt`):

| | Sweep | `RUVIII_C_Varying` | RUV worker | Normalisation |
|---|---|---|---|---|
| Metabolomics | 114 s | 467 s | 602 s | 789 s |
| Lipidomics | 204 s | 830 s | 1,065 s | 1,343 s |

- **The third-party correction is 78% of RUV and about 60% of normalisation.** It is `RUVIIIC::RUVIII_C_Varying`, compiled code that proteomics also calls. Nothing here changes it.
- **The sweep is 19% of RUV.** Running the ANOVA only once could save at most part of that. Proteomics does not do this either, so it is outside this ticket.
- **Lipid is 1.8 times metab in both parts.** The same factor in both suggests the host was busier during `084d`; this was not established.

## With missing values

One NIPALS imputation at 200 × 18,000 with 5% missing did not finish inside a 1,200 s cap, so at that size:

**Corrected 2026-11-05 (KIT-086):** the two estimates below came from uniform missingness, which the corpus doesn't use. On the declared clustered pattern, old imputation measured about 5 minutes per assay, and KIT-087 then removed that cost. RUV fits the deadline at 603-617 s. The lines stay as written.

- **Before:** 40 imputations over 2 assays, at least 13 hours.
- **Now:** 2 imputations, at least 40 minutes.

These are lower bounds from one simulated, capped run on a shared host. No journey has confirmed them.

**The worker deadline remains a limit.** At least 40 minutes of imputation exceeds the application's 1,800 s deadline per worker step, so this RUV step would still be stopped at this size with this much missing data. The change makes it possible to fix; it does not make it fit.

## Same standard as proteomics: what was not brought across

- **The viability rule (`ruvCompleteControlRows`).** Proteomics skips a percentage whose elbow k exceeds its complete control rows (KIT-005, a user decision after a correction returned no finite values). The omics sweeps do not. It can change which percentage is chosen, so it is a scientific change and needs its own decision.

## Divergences

- **Screen journeys in place of large ones,** as above, and told to the user before running.
- **The correction pass was used.** The first launches ran no member: my binding edit failed, and `run.sh` carried on. Recorded in `power-state.txt`; the repair is 1d8c64f6.
- **The ticket's shape differs from proteomics':** the sweep passes a prepared object through `ruvCancor` rather than calling a from-input function, to keep the S4 generic's signature and the tests' pinned call count of one per percentage.
