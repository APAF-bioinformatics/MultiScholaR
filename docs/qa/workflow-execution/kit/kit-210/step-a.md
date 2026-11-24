# KIT-210 step A: what RUV-III-C does with a missing value (A1 to A6)

Read at the step-A commit e8a43097 (v0.9-staging c03cb02a plus this ticket's records; `git diff --stat c03cb02a
53713252 -- R/ src/ inst/ packaging/` is empty, so the reading also holds at the staging head of 2026-11-24).
Upstream source: RUVIIIC 1.0.19, `cran/RUVIIIC` at c534151f (the lock's commit), read from the archive whose
member hashes A1 re-computed on a runner. Marks: **shown** (read at the file and line cited, measured in the run
cited, or an algebra step written out here), **inferred** (derived, not observed), **conjecture** (step C tests it).

## A6. The answer to the owner: "can RUV-III-C safely run on zeros?"

In five lines, each marked:

1. **Shown (source, and the installed library on a runner):** the zeros are RUVIIIC's own. Every form of the
   library sets each missing value to 0, always, with no option (S2), and those zeros act only in the step that
   finds the k unwanted directions (S6, A-1). The installed library does exactly this (A2, run 45351569441). The
   package passes missing values as NA, as the library asks, and never writes a zero itself (A4). The library offers
   no other handling of missing values (A-I).
2. **Shown (the library's own words and test):** its documentation says zeros "appear to be much larger than the
   discrepancies due to other (correctable) technical factors" so that RUV-III "will attempt to correct for the
   larger (uncorrectable) discrepancy" (`R/RUVIII_C.R:10`); it avoids this for the feature being corrected only. Its
   test of the zero-filling uses standard-normal data, where a zero is the mean (`test-compare_RUVIII_Varying.R:37`),
   not log-intensities, where a zero sits about 20 log2 units from the value it replaces.
3. **Shown (algebra, A-3, A-4):** the zeros distort the matrix the directions come from in three parts: each run's
   zero mass, the zero mass each pair of runs shares, and a smaller cross part. Even if every run had the same number
   of zeros, only the first part cancels; the directions stay put only under a pattern random missingness does not
   produce. **Measured (KIT-185, VAL-140):** at random missingness from 0.5 % (P2) or 2 % (most shapes) the
   zero-filling leaves 0.09 to 0.40 of the unwanted variance.
4. **Inferred (first-order perturbation, A-5):** the error grows with the squared log level and the missing rate,
   and falls with the strength of the unwanted factors; the shared-zeros part grows with the number of runs times
   the rate. No safe threshold follows from the algebra.
5. **Conjecture (step C tests it):** that a quantity computed from a user's own data (kappa-hat) can mark the data
   sets where the zero-filled form is safe; that mean filling (C6) or NIPALS filling (C7) of the direction step
   removes most of the damage.

**Plain answer:** "On log-scale data, not in general. The zeros are the library's, not ours; they sit about 20
log2 units from the values they stand for, and even when spread evenly they turn the unwanted-variation directions
at the missing rates real data have. Whether a check on the user's own data can mark the safe cases is what the
simulation (step C) tests." This confirms ticket 2.7; nothing in it is overturned. Not read: Poulos et al. 2020
(A3), so whether its authors state conditions for the zero-filling is not known.

## Outcome of step A

**A-I.** GA1 passes (A1: the CRAN archive and the GitHub mirror at c534151f are identical file by file, 32 of 32,
and equal to ticket 2.2; the installed package is Version 1.0.19, RemoteSha c534151f). GA2 passes (A2: every
deciding part (a), (b), (c1), (c2), (c3) passes, with its negative controls). GA3 is not applied (A3 not read).
The library has no missing-value handling besides the zero-filling (S1, S2, S8; the installed R functions equal
the source's, `functions.tsv`). A-II: no such form in the source. A-III: GA2 passed. A-IV: the package still
passes NA and fills nothing before the call (A4).

## A1. The pinned source, fetched and hashed on a runner

Run 45351569441 (request `k210-step-a-e8a43097`, package commit e8a43097, reference BLAS, R 4.6.1). Files
under `step-a/`: `tarballs.tsv`, `source-hashes.tsv`, `md5-check.tsv`, `installed.tsv`, `functions.tsv`,
`summary.txt`, `console.txt`, `session.txt`, `meta-45351569441.txt`.

| Check | Observed (shown) |
|---|---|
| CRAN archive | `src/contrib/RUVIIIC_1.0.19.tar.gz` returns 404 (the package is no longer current on CRAN); the archive copy `src/contrib/Archive/RUVIIIC/RUVIIIC_1.0.19.tar.gz` was used, SHA-256 c6556d081795536c10c3180766826f23645e6748cc639d15eff3b56ac93da8df (399,491 bytes) |
| GitHub mirror | `cran/RUVIIIC` at c534151f, SHA-256 846cc1442ad90711efade03c9e2f6e90c963753023860a91efe9ef63958d4009, equal to the scratchpad tarball the ticket read |
| File by file | 32 files in each archive, the same 32 names, every SHA-256 equal (`source-hashes.tsv`); each archive's own MD5 list matches its 31 listed files |
| Ticket 2.2 | the six recorded hashes match in both archives, including `src/RUVIIIC.cpp` 57615f70...48c0e1797 (the r1.1 correction) |
| Installed | Version 1.0.19; RemoteType github, cran/RUVIIIC, RemoteSha c534151f; built on the runner from that source (`packaging/r/install-library.R:33`); `libs/RUVIIIC.so` SHA-256 13510d12...b2d4 (a build on that runner, not a reference value) |
| Installed R layer | `RUVIII_C`, `RUVIIIC_R`, `RUVIII_C_Varying`, `RUVIIIC_Varying_R`, `RUVIII_C_Varying_residual_dimension`: formals and bodies identical to the source's (reported) |

GA1: **pass**.

## A2. The installed library zero-fills, and only in the direction step

Same run. A 12 x 60 log-scale matrix (levels N(20, 2^2), two unwanted factors orthogonal to three replicate
sets of four runs, a set effect on 22 features, noise sd 0.2), target f01 missing at r03 (11 rows, residual
dimension 8, so the eigen branch), 10 complete potential controls, six missing values in features that are
neither; k = 2; one OpenMP thread. `step-a/a2.tsv`:

| Part | Rule | Observed | Pass |
|---|---|---|---|
| (a) the NA of f20 at r05 (neither target nor control; the target is observed at r05) replaced by 0 | output `identical()` | identical | yes |
| (b) the same NA replaced by 20 (negative control of (a)) | differs by > 1e-6 | 5.5e-02 | yes |
| (c1) y - W alpha(target) against the returned corrected values | <= 1e-10 | 0 (and NA on the missing row) | yes |
| (c2) W, the controls' alpha and W alpha(target) from U = top-2 eigenvectors of P Sel (Y0 Y0') Sel' P computed in R from the zero-filled matrix, basis aligned | <= 1e-8 | 4.0e-13, 3.5e-12, 4.3e-13; R'R - I 1.3e-12; relative eigengap 0.50 | yes |
| (c2-neg, reported) the same reproduction from the matrix filled with 20 | > 1e-6 | 0.25 | yes |
| (c3) 0: W, alpha, residual dimension, output identical; 20: U and output change | 0 identical; 20 > 1e-6 | identical; U projector 0.11, output 0.055 | yes |
| (d, reported) `RUVIII_C` (fixed controls): 0 identical, 20 differs | as (c3) | identical; 0.055 | yes |

Read: a missing value of a feature that is neither target nor control is, to the compiled code, exactly a 0
(a, c3); it changes the result through the directions (b, c3: U moves), and the result is reproduced exactly by
RUV-III on the zero-filled cross-product (c1, c2), not by a different fill (c2-neg). GA2: **pass**.

## A3. The paper (Poulos et al. 2020, doi:10.1038/s41467-020-17641-3)

**Not read.** The route of reading is the coordinator's E4, which is not decided; the brief forbids calling an
outside service from this machine, and the runner job may download only the pinned RUVIIIC source. So GA3 is not
applied: it neither passes nor fails. **What therefore stays unverified:** whether the method's authors state
conditions under which the zero-filling of the other columns is acceptable (for example a missing rate, or that
the data are centred), and whether they describe any alternative handling of missing values in the article that
the library does not ship. A-I does not depend on it: the library's own options are settled by its source (A1,
S1, S2, S8); what the paper could add is a stated assumption for the decision table (GA3's consequence).
The library itself cites the paper (`DESCRIPTION:12-13`) and states the rationale in its documentation
(`R/RUVIII_C.R:8-12`), quoted under A4/S3 below.

## A4. The library and the package at the step-A commit (file:line)

### The library (RUVIIIC 1.0.19 at c534151f)

| # | Shown | Where |
|---|---|---|
| S1 | Exports: `RUVIII_C`, `RUVIII_C_Varying`, `RUVIII_C_Varying_residual_dimension`, `omp_get_num_threads`, `omp_set_num_threads`; nothing else. Each correction function has a `version = "CPP"` (default) and `"R"` form | `NAMESPACE:3-7`; `R/RUVIII_C.R:52-78`; `R/RUVIII_C_Varying.R:48-74` |
| S2 | Every form copies the whole matrix, sets each missing value to 0, and forms the runs-by-runs cross-product once | `src/RUVIIIC_Varying.cpp:80` (copy), `:83` (`std::isnan(v) ? 0.0 : v`), `:85` (Y0 Y0'); `src/RUVIIIC.cpp:88`, `:91`, `:93`; `R/RUVIII_C_Varying.R:110-111`, `:127`; `R/RUVIII_C.R:126-127`, `:143`. No argument reaches these lines (signatures `RUVIIIC_Varying.cpp:25`, `R/RUVIII_C_Varying.R:48`, `R/RUVIII_C.R:52`) |
| S3 | The documentation states the harm and the rationale: zero values "will appear to be much larger than the discrepancies due to other (correctable) technical factors. RUV-III will attempt to correct for the larger (uncorrectable) discrepancy, and ignore the correctable technical factors" (`:10`); "other columns besides X will still have missing values. These values are replaced with zero in order to apply RUV-III. No additional transformation is applied to the input data matrix" (`:12`) | `R/RUVIII_C.R:8-12` |
| S4 | Input: Y "should contain missing (NA) values, rather than zeros" | `R/RUVIII_C_Varying.R:13`; `R/RUVIII_C.R:17` |
| S5 | Per target: the rows where it is observed (`.cpp:122-133`); potential controls observed on every one of those rows (`:138-150`); replicate sets with a row; residual dimension rows minus sets (`:151-166`) | `src/RUVIIIC_Varying.cpp` |
| S6 | The target's rows are copied from the zero-filled matrix (`:127`); directions = top-k eigenvectors of P Sel (Y0 Y0') Sel' P (`:241`, Spectra largest magnitude `:248`); alpha = U' Y0[S, ] for every column (`:270`); the chosen controls' alpha (`:278`) and values (`:279`); W = Yc ac' (ac ac')^-1 (`:281`); corrected = y - W alpha(target) (`:284`). Only the target's and the chosen controls' columns of alpha are used, and both are observed on every row of S (S5), so **no zero reaches the correction except through U** | `src/RUVIIIC_Varying.cpp` |
| S7 | Fewer than k residual dimensions or controls: NA (`:168-171`); exactly k: a GLS formula on the controls only, no zero (`:173-198`) | `src/RUVIIIC_Varying.cpp` |
| S8 | `RUVIII_C` stops on a missing control (`src/RUVIIIC.cpp:67-73`; `R/RUVIII_C.R:101-104`) and zero-fills the other columns too (S2) | as cited |
| S9 | The compiled functions check only `k >= nRows` and `k > controls` (`RUVIIIC_Varying.cpp:44-51`); the R forms refuse `k <= 0` (`R/RUVIII_C_Varying.R:81-84`) | as cited |
| S10 | Upstream's test of the zero-filling uses `rnorm(20*20)` data (mean 0) and compares with `ruv::RUVIII` on the data with the other column's NA set to 0 | `tests/testthat/test-compare_RUVIII_Varying.R:37`, `:58` |

No form, option or argument handles a missing value other than by this zero-filling (S1, S2, S8; a search of
`R/`, `src/` and `DESCRIPTION` for "averag", "impute" and "missing" finds only the zero-fill lines, the row
selection and the documentation). The DESCRIPTION's "applying RUV-III to complete subsets of the data"
(`DESCRIPTION:12`) describes the row selection per target (S5), not the columns: the other columns are zero-filled.

### The package: what it passes (re-read; matches ticket 2.3)

| Route | Call | What a missing value is when passed | Shown at |
|---|---|---|---|
| R1 proteomics legacy | `diaNormalizationOperation` `R/utils_workflow_dia_products.R:490` -> S4 method `R/func_prot_s4_norm_methods.R:248-254` | NA: Y is `protein_quant_table` transposed (`:228`); the normaliser step before it sets every non-finite value to NA (`:69`, `:86`). The RUV method itself does not re-check (inferred: a value that is non-finite but not NA could only arrive if the RUV step were run without the normaliser step) | k < 1 refused `:212` -> `R/utils_workflow_dia_norm.R:34-45`; controls guard `:231-235`; levels restored `:262-263` |
| R2 proteomics declared | `R/utils_workflow_dia_products.R:479-480` | NA, as R1 | DA uses W from `ruv::RUVIII` on complete proteins (`R/func_prot_norm_optimization_helpers.R:1352-1357`) |
| R3 lipid / metab legacy | `R/utils_workflow_omics_declarations.R:355-356` -> `R/func_lipid_s4_ruv_methods.R:691-697` ("Use data potentially containing NAs"); metabolomics identical lines (`diff` of lines 470-720 with the class names normalised is empty) | NA: `:549` sets non-finite to NA; "NO IMPUTATION" `:566`; `:572` warns that NAs are present | `:678-682` controls guard; k >= 0 accepted `:481` |
| R4 lipid / metab declared | `R/func_omics_ruv_factors.R:361-362`, levels restored `:371` | NA | W from `ruv::RUVIII` on complete features `:137-176`; k = 0 returns "none" `:147` |

`grep` over `R/` finds no other `RUVIII_C_Varying` call than these and the peptide method
`R/func_pept_s4_norm_methods.R:883` (no caller, as ticket 2.3). Import turns zeros into NA for lipid/metab
(`R/utils_workflow_omics_import.R:369`). **The package never writes a zero and does not fill a value before the
call (A-IV does not hold).**

### The route inconsistency: k chosen on NIPALS-imputed data, correction on zero-filled data (shown)

- Proteomics, Automatic: `resolveProtRuvParameters` calls `findBestNegCtrlPercentage` (`R/utils_workflow_dia_norm.R:68-76`),
  whose input is `prepareRuvCancorInput` (`R/func_prot_norm_optimization_helpers.R:793`), which NIPALS-imputes when
  the matrix has a missing value (`R/func_prot_s4_norm_methods.R:169-181`); `best_k` becomes k
  (`R/utils_workflow_dia_norm.R:77-78`); the correction is RUV-III-C on the matrix with its NAs, zero-filled inside
  the library (`R/func_prot_s4_norm_methods.R:248-254`; S2).
- Lipidomics / metabolomics: the sweep's canonical-correlation object is `prepareOmicsRuvCancorObject`
  (`R/func_lipid_norm_ruv_optimization_helpers.R:87-89`; `R/func_metab_norm.R:412-414`), which NIPALS-imputes
  (`R/func_metab_norm.R:293-318`, the call `:308`); the `ruvCancor` S4 method imputes the same way
  (`R/func_lipid_s4_ruv_methods.R:150-165`); the correction is RUV-III-C (`:691-697`).
- **Correction to ticket 2.3 and section 5 (C7).** The number of NIPALS components is not `inst/config/config.ini:87`:
  that line is in the `[ruvCovariates]` section (the limpa route's covariates), not the legacy RUV-III-C route.
  What the legacy routes use: proteomics sweep 5, the default argument of `findBestNegCtrlPercentage`
  (`R/func_prot_norm_optimization_helpers.R:679`; `resolveProtRuvParameters` does not pass it); the proteomics
  `ruvCancor` S4 method 2 by default (`R/func_prot_s4_norm_methods.R:135`); the lipid/metab sweep the object's stored
  `ruvCancor` argument, else 2 (`R/func_metab_norm.R:298`); the lipid/metab `ruvCancor` S4 method 2 by default
  (`R/func_lipid_s4_ruv_methods.R:28`). C7's `ncomp = 5` therefore matches the proteomics Automatic sweep, not the
  lipid/metab default (2). This does not change A-I; it is a fact for C7's definition in the declaration.

## A5. The algebra, self-contained

Notation, one target j. S: the m runs where j is observed. Y*: the m-by-n matrix of the values that would have been
observed on S (all features). B: the missingness indicator (B_ij = 1 when feature j is missing at run i).
Y0 = Y* with every missing entry set to 0, which is what `.cpp:83` makes, restricted to S by Sel (`.cpp:127`,
`:241`). y*_ij = mu_j + d_ij (feature level plus deviation). M_S: the replicate matrix on S with the sets that have
a run (`.cpp:151-166`). P = I - M_S (M_S' M_S)^-1 M_S'. U: the k leading eigenvectors of P Y0 Y0' P (`.cpp:241-268`).
Yc: the chosen controls on S (observed everywhere on S). A_c = U' Yc. W = Yc A_c' (A_c A_c')^-1 (`:281`).
Corrected y^ = y_S - W U' y_S (`:270`, `:284`).

Conditions: (i) every run of S lies in exactly one replicate set, so 1 = M_S 1 and P 1 = 0 (M is
`model.matrix(~ grouping - 1)`; dropping rows and empty sets keeps one 1 per row); (ii) U is the eigenvector set of
the top k eigenvalues with a gap, lambda_k > lambda_{k+1}; (iii) A_c has full rank k; (iv) a perturbation bound is
informative only when the perturbation is small against the eigengap.

1. **A-1 (exact, from S6).** y^ depends on the zero-filled entries only through U: y_S and Yc are observed on S,
   W is a function of Yc and U, alpha(target) = U' y_S.
2. **A-2 (exact; needs (iii)).** U' W = A_c A_c' (A_c A_c')^-1 = I_k, so (W U')^2 = W (U' W) U' = W U': y^ =
   (I - W U') y_S is an oblique projection that removes the part of y_S in span(W) = span(Yc Yc' U) along the null
   space of U'. If U is rotated away from the unwanted space, the wrong part is removed and unwanted variation stays.
3. **A-3 (exact; needs (i)).** (Y0 Y0')_{ii'} = sum_j y*_ij y*_i'j (1 - B_ij)(1 - B_i'j). Expanding the product:
   Y0 Y0' = Y* Y*' - A - A' + Q, with A_{ii'} = sum over j missing at i' of y*_ij y*_i'j, and Q_{ii'} = sum over j
   missing at both i and i' of y*_ij y*_i'j (Q = E E', E holding Y* at the missing entries; c_i = Q_ii is run i's
   zero mass). Split A_{ii'} = v_i' + R_{ii'} with v_i' = sum over j missing at i' of mu_j y*_i'j (depends on i'
   only: A's level part is the matrix 1 v') and R_{ii'} = sum over j missing at i' of d_ij y*_i'j. By (i),
   P (1 v') P = (P 1) v' P = 0 and P (v 1') P = 0. Hence **P Y0 Y0' P = P Y* Y*' P + P (Q - R - R') P, exactly.**
4. **A-4 (exact for what it states).** If every run carries the same zero mass, c_i = c, the diagonal of Q adds
   P (c I) P = c P, which on range(P) shifts every eigenvalue by c and moves no eigenvector. **Only the diagonal part
   is removed.** U is unchanged exactly when P (Q - R - R') P is a multiple of P on range(P). A sufficient condition
   is Q = a I + b 1 1' (equal mass on every run and equal shared mass on every pair) together with P (R + R') P = 0;
   then P Q P = a P. Equal per-run mass does not give it: the off-diagonal (co-missing) entries of Q and the
   projected cross part R remain. Random masks give neither (inferred: under MCAR at rate r the co-missing counts
   vary between pairs with sd about r sqrt(n)).
5. **A-5 (inferred, first order; needs (ii), (iv)).** Davis-Kahan (or the Yu-Wang-Samworth form): sin of the largest
   angle between U and the unwanted space <= C ||N||_op / delta, N = the non-isotropic part of P (Q - R - R') P,
   delta = lambda_k - lambda_{k+1} of P Y* Y*' P (about n tau^2 ||P z||^2). Orders under MCAR (inferred): the
   per-run part acts with size about mu^2 sqrt(n r (1 - r)), the co-missing part about mu^2 r sqrt(n) sqrt(m); R is
   smaller by about sd(d)/mu unless the missingness selects on d (left-censoring does). Hence kappa = (mu^2 / tau^2)
   sqrt(r (1 - r) / n) sqrt(1 + r m) / ||P z||^2, with ||P z||^2 = m - 1 when the replicate sets are the groups and
   the factors are orthogonal to them (P5, sets = units: about m - 8). At KIT-185's harmful rates (iv) fails, so this
   is a description of direction, not a law. **Only steps 3 and 4 are exact; nothing here shows a safe threshold.**
6. **A-6 (exact on complete data; consequence inferred).** Subtracting a constant per feature (c_j) before the
   call: Y -> Y - 1 c'. By (i) P (Y - 1 c')(Y - 1 c')' P = P Y Y' P, so U is unchanged; U lies in range(P) and P 1 = 0, so
   U' 1 = 0 and alpha = U' (Y - 1 c') = U' Y; W changes by
   1 w' (w' = c_c' A_c' (A_c A_c')^-1); y^ changes by a constant per feature, which `restoreRuvFeatureLevels`
   (`R/func_prot_norm_optimization_helpers.R:1133-1140`) removes on R1/R2 (and `restoreOmicsRuvAssayLevels`, `R/func_omics_ruv_factors.R:180`, on R4); on R3 (no restoration) the
   result differs from S's by a constant per feature.
   With missing values, the fill after centring stands for the feature's observed mean (candidate C6, a remedy; it
   does not make the zero-filled form safe).
7. **A-7 (inferred).** The value-dependent false-discovery rise that S and C share in KIT-185 comes from testing
   incomplete features with the legacy test's degrees of freedom (VAL-53 (1)), not from the directions; A-3 to A-5
   do not explain it.
8. **A-8 (definition).** kappa-hat as ticket 2.4 A-8 (the co-missing term included).

## Not checked in step A

- Poulos et al. 2020 (A3, above).
- The library's `version = "R"` forms on the runner (their zero-filling is shown from source, S2; A2 ran the
  compiled forms, which are the package's calls). The installed R wrapper functions were compared with the source
  (`functions.tsv`).
- Whether k = 0 can reach the legacy lipid/metab method in practice (ticket 2.9; unchanged).
