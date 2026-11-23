# KIT-183 feasibility gates

## F1: limma 3.68.4 equals 3.68.5 as the auditor read it (read, no run)

Observed (2026-11-23, by the implementer, on the laptop; no R run):

- Tarball `limma_3.68.4.tar.gz` from `https://bioconductor.org/packages/3.23/bioc/src/contrib/Archive/limma/`,
  SHA-256 `7bfc3cd72754dad6e5edb2ef7c81f5ae85837a34375269e69738801d6e99a210`; its `DESCRIPTION` reads
  `Version: 3.68.4`, `git_last_commit: 7c11dae` (the lock's `GitLastCommit`, `packaging/multischolar.lock.json`).
- Tarball `limma_3.68.5.tar.gz` from `https://bioconductor.org/packages/3.23/bioc/src/contrib/`, SHA-256
  `40b1d2059a76ad98db4bf80139760acf8e8bae92a8265aaaf16ac199831533f4`.
- `R/geneset-camera.R`, `R/geneset-cameraPR.R`, `R/zscore.R`, `R/squeezeVar.R`, `R/ebayes.R`: byte-identical between
  the two (`cmp`); SHA-256 values below.
- Read in 3.68.4's `R/geneset-camera.R` (header "Last modified 4 Mar 2025"), each fact of ticket 2.2:
  - degrees of freedom `:72-80`: preset correlation `G-2L` (ranks `Inf`); estimated `min(df.residual, G-2L)`;
    cameraPR `G-2L` (`geneset-cameraPR.R:47-51`);
  - `trend.var` `:146-147`: `A <- rowMeans(y)` passed as `covariate` to `squeezeVar`;
  - the statistic `:148-155`: `modt <- unscaledt / sqrt(sv$var.post)`, `zscoreT(modt, df = df.total, approx = TRUE,
    method = "hill")`, `df.total <- min(df.residual + sv$df.prior, G*df.residual)`;
  - `varStat <- var(Stat)` over all G features `:159-160`; `delta <- G/m2*(meanStatInSet-meanStat)` (the set's mean
    minus the other features' mean), `varStatPooled`, `two.sample.t` `:193-196`; `Direction` from the sign `:205-209`;
  - the reported correlation `:186-192`: `tab[i,2] <- correlation` is written before `if(!allow.neg.cor) vif <-
    max(1,vif)`;
  - `directional = FALSE` forces ranks `:65-70`;
  - a character `contrast` is matched to the design's column names, an integer below `p` moves that column last
    `:101-114`.
- `eBayes(trend = TRUE)` (`R/ebayes.R`): covariate `fit$Amean`, `squeezeVar(sigma^2, df.residual, covariate = ...)`,
  `df.total <- pmin(df.residual + df.prior, sum(df.residual))`: for complete data the same moderated t and `df.total`
  camera computes internally (used for `effect` and for cameraPR's statistic).

Gate result: **passes**. No difference in any fact 5.1 relies on.

Full SHA-256 of the files read (3.68.4 = 3.68.5):
- `R/geneset-camera.R` `4c4c0b437f05f0fab3ab38ae2d865f7e20fd7364067917b71e44f790883dc23f`
- `R/geneset-cameraPR.R` `ecfae11557f45b4467453c5f2f0dd27a5e02ee7918390e2478057de474a75c9c`
