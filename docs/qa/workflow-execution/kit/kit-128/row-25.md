# KIT-128: parked row 25 recorded (decision D5; no code)

Written 2026-11-17 on branch `kit-128`. Nothing run for this record; the reading is W0's section 6 (`w0-design.md`).

**Row.** DIA_LIMPA on data with no missing values: `limpa::dpcDE` can stop with "No residual degrees of freedom in
linear model fits"; 2 of 5 runs in the KIT-124 S2 pilot, cell S5.

**State.** Upstream, in limpa (1.4.0 installed; the same lines in the current release 1.4.2: `R/dpcDE.R:7-8`,
`R/voomaLmFitWithImputation.R:136-139` and :234-237). Our call (`R/func_prot_da_model_methods.R:391-392`) passes the
documented arguments.

**Mechanism (from the code).** With no missing value, `dpcDE` passes a predictor that is the log of the DPC standard
errors. When those are one constant (complete data, the same precursor count for every protein, variances squeezed
to a common value), the predictor's row means are collinear with the intercept of the variance-trend fit; `lm.fit`
drops the column and returns `NA` for its coefficient; every weight is `NA`; `lmFit` treats every observation as
missing; `limma::eBayes` (our converter, `R/func_prot_limpa_da_helpers.R:52`) stops with the message. Whether it
happens depends on floating-point noise in the standard errors (`lm.fit`'s tolerance 1e-7), hence 2 of 5 runs.

**Decision D5.** Recorded only. To be reported to the limpa authors from the user's account (the user's call). No
guard on our side unless it shows on real data, where varying precursor counts make the standard errors vary.
