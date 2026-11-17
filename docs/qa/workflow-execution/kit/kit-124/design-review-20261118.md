# KIT-124: which RUV designs to offer (literature and our simulations, 2026-11-18)

Written for the user's gate, on the user's questions of 2026-11-18: can every design be an option with a warning; why does a literature-faithful implementation fail; are pseudo-replicates acceptable; is "all features as controls" correct, or should controls be features found null by a first test. Nothing here changes code or a default.

## 1. What the method papers say

| Question | Source | What it says |
|---|---|---|
| Controls for replicate-based RUV (RUVs) | RUVSeq vignette, RUVs section | "Although in principle one still needs control genes ... RUVs is robust to that choice and ... using all the genes works well in practice." Its code: `RUVs(set, genes, k=1, differences)`, `genes` = all genes, `differences <- makeGroups(x)` with x the treatment groups. |
| Controls for control-only RUV (RUVg) | RUVSeq vignette, RUVg section | A priori controls (spike-ins) or "in-silico empirical" controls: "least significantly DE genes based on a first-pass DE analysis"; its code keeps "all but the top 5,000 genes" by edgeR p-value (most genes kept, the clearly changed ones dropped). |
| Test on corrected data or factors in the model | RUVSeq vignette | Normalised counts "should be used only for exploration"; "It is important that subsequent DE analysis be done on the original counts ... as removing the unwanted factors from the counts can also remove part of a factor of interest"; design `~ x + W_1`. |
| The same, in general | limma, `removeBatchEffect` help | "This function is not intended to be used prior to linear modelling. For linear modelling, it is better to include the batch factors in the linear model." |
| The same, measured | Nygaard, Rødland, Hovig 2016, Biostatistics 17:29 | Two-step correction then testing "may lead to exaggerated confidence in downstream analyses", worst when groups and batches are unbalanced. |
| RUV-III: controls | Molania et al. 2019, NAR 47:6073 | "the question will not be 'Can all genes be suitable negative controls?' but 'Does using all genes as negative controls help?' We will see examples where doing so does help, and others where it does not. Selecting suitable negative controls is the major challenge in using RUV-III." |
| RUV-III: replicates | Molania et al. 2019 | Makes "vital use of technical replicates and suitable control genes"; "an approach using pseudo-replicates when technical replicates are not available". No clear rule for k. |
| RUV-III: the test | Molania et al. 2019; Jiang, Gagnon-Bartsch, Speed (arXiv 2310.01918) | Both run limma on the RUV-III-adjusted matrix. Jiang et al.: "the asymptotic theory developed in this paper does not address the validity of the p-values"; "the choice of negative controls and the specification of what constitutes a replicate plays a critical role". |
| RUV-III-C (proteomics) | Poulos et al. 2020, Nat Commun 11:3793 | True technical replicates across instruments; negative controls known a priori (yeast peptides, scaled); k = 20 on about 1,500 runs; no differential test on the corrected data in the paper. "Key to success ... technical replicates ... and also a suitably large number of negative controls." |
| Pseudo-replicates by known biology | Salim et al. 2022, NAR (RUV-III-NB) | "cells with the same level of the biological factor of interest can be declared to be a pseudo-replicate set"; a heterogeneous set "may risk removing some of the biological signal of interest". Recommends housekeeping genes as controls; performance "relatively robust to ... the choice of negative control gene sets"; the one set chosen as non-DE from the data was the only one affected by overestimating K. |
| Metabolomics | Kim et al. 2021, Nat Commun 12:4992 (hRUV) | Data-driven non-DE controls against all metabolites as controls: "no significant differences between the two approaches"; sample replicates "were more effective for batch correction". |

Reading: (a) groups as pseudo-replicates are a literature convention (RUVSeq's own RUVs example; RUV-III; RUV-III-NB). (b) For replicate-based RUV, all features as controls is an accepted, documented fallback when no a-priori controls exist; a-priori known controls are the first choice. (c) First-pass "null by a test" controls are the convention for control-only RUV (RUVg), where the clearly changed features are dropped and most are kept. (d) The RNA-seq and limma convention is factors in the model; the RUV-III papers test the corrected matrix, and their own theory paper says the p-values' validity is not addressed.

## 2. Where today's MultiScholaR designs sit

The computations match the packages (KIT-123 W1 reference tests). The usage differs from the papers in three choices, in every omic:

1. **Controls:** the least significant 5 to 30% of features by a first test, chosen by a grid (`R/func_trx_ruv.R:5,18`; the proteomics sweep). The vignette's RUVs passes all genes; its RUVg keeps all but the most significant.
2. **k:** chosen by the grid (transcriptomics mean 2.9 where the truth is 1); the vignette uses 1.
3. **The test:** classic proteomics, lipidomics and metabolomics test the corrected matrix (the RUV-III papers' usage); transcriptomics and DIA_LIMPA put W in the model (the vignette's).

Groups as replicate sets (pseudo-replicates) are today's replicate structure and are also the replicate structure of every design that holds below.

## 3. What our simulations measured (core cells; nominal 0.05; worst cell with true effects)

Mean false discovery proportion, worst core cell (and in C5, factor tied to groups at rho 0.3); power in the centre cell.

| Design | Classic prot | DIA_LIMPA | Transcriptomics | Lipid/metab |
|---|---|---|---|---|
| No RUV | 0.052; power 0.49 | 0.191; 0.56 | 0.044; 0.005 | 0.059; 0.003 |
| True factor known (ceiling) | 0.052; 0.91 | 0.055; 0.91 | 0.055; 0.77 | 0.061; 0.50 |
| Today's design (arm 3) | 0.385; 0.92 | 0.466; 0.91 | 0.299; 0.73 | 0.556; 0.42 |
| Today's controls, W in the model (arm 4) | 0.322; 0.89 | 0.435; 0.90 | (arm 3) | 0.162; 0.17 |
| First-pass nulls: drop the most significant 25%, W in the model (arm 5) | 0.078; 0.91 | 0.149; 0.92 | 0.064; 0.60 | 0.239; 0.50 |
| All features, groups as replicates, W in the model (arm 7) | 0.066; 0.91 | 0.055; 0.91 | 0.055; 0.77 | 0.072; 0.43 |
| Known true nulls, W in the model (arm 8) | 0.053; 0.91 | 0.056; 0.91 | 0.056; 0.76 | 0.068; 0.37 (standards) |
| Known nulls, test the corrected matrix (arm 10) | 0.120; 0.93 | 0.152; 0.93 | 0.072; 0.79 | 0.267; 0.57 |
| Supplied set with 20% true effects in it (arm 9) | 0.119 | 0.089 | 0.076 | 0.371 |

- **Arm 8 against arm 10** isolates the test form (same controls, same k): testing the corrected matrix roughly doubles to quadruples the false discovery proportion in three of four omics.
- **Arm 4 against arm 7** isolates the control rule (both W in the model, both pseudo-replicates): the small least-significant subset is what fails, mostly when the factor is tied to the groups.
- **Arm 8** shows the user's premise is right: true nulls as controls are the best design. **Arm 5** shows the cost of finding them from the same data by a first test: acceptable in classic proteomics and transcriptomics except when the factor is tied to the groups (and with a power loss in transcriptomics), poor in DIA_LIMPA and lipid/metab.

Stress cells finished so far (transcriptomics, lipid/metab; scored 2026-11-18, partial):

| Cell | Transcriptomics: arm 7 / arm 5 / today | Lipid/metab: arm 7 / arm 5 / today |
|---|---|---|
| 40% of features changed (S2) | 0.035 (power 0.83) / 0.035 (0.09) / 0.036 | 0.055 / 0.066 / 0.318 |
| Factor strongly tied to groups, rho 0.6 (S3) | 0.052 / 0.062 (power 0.00) / 0.724 | 0.130 / 0.537 / 0.633 |
| 80% of effects in one direction (S9) | 0.074 (= ceiling 0.074) / 0.058 / 0.074 | 0.059 / 0.096 / 0.356 |
| Loadings tied to abundance (S10) | 0.048 / 0.049 / 0.047 | 0.049 / 0.095 / 0.430 |
| Three factors (S4): arm 7 with k estimated / k = 1 | 0.049 (power 0.76) / power 0.00 | 0.037 / n.a. |
| 3 per group (S1) | 0.046 / 0.044 / 0.046 | 0.101 (ceiling 0.126) / 0.070 / 0.546 |

All features as controls held where it was expected to be weakest (many changed features, one-directional effects, abundance-tied loadings) in transcriptomics; in lipid/metab it is liberal under strong confounding (0.130 against a ceiling of 0.072) and at 3 per group, where the ceiling is liberal too. k = 1 fails with three factors: k must be estimated or chosen, not fixed. Classic proteomics and DIA_LIMPA stress cells are still running.

## 4. Proposal for the option set (for the user's decision)

Three tiers, each a tested arm, each recorded in the analysis record and named in the report's methods text; the measured figure is shown beside the option, not only in a tooltip.

| Tier | Option | Label shown | Omic |
|---|---|---|---|
| Recommended | Pseudo-replicates (groups), all features as controls, factors in the model (arm 7) | "Recommended" | all four; the default where the user accepts it |
| Valid with conditions | Technical replicates as the replicate structure (arms 15/16) | "When technical replicates are declared" | all |
| Valid with conditions | Supplied controls (housekeeping list, internal standards), factors in the model (arm 8) | "When you have a control list; a list contaminated with changing features inflates false discoveries" | all; standards for lipid/metab |
| Valid with conditions | First-pass null controls (drop the most significant 25%), factors in the model (arm 5) | "Not when the unwanted factor follows the groups; loses power in transcriptomics" | classic proteomics, transcriptomics |
| Legacy, warned | Today's design (arm 3) | "Legacy (results before version X); measured false discovery proportion 0.30 to 0.56 at a nominal 0.05" | all |
| None | No RUV | "Calibrated except DIA_LIMPA; loses most power when unwanted variation is present" | all |

The corrected matrix is still computed and used for PCA, plots, export and integration under every RUV option: "factors in the model" changes only what the test reads.

Open for the stress cells still running: conditions of use for classic proteomics and DIA_LIMPA; whether k is estimated (`sva::num.sv`, not a dependency today) or user-chosen.

## 5. Choosing k: two rules, offered where each is defined (user's direction, 2026-11-18)

| | Estimated count (`sva::num.sv`, method "be") | Separation rule (today's `findBestK`, `trxFindBestK`) |
|---|---|---|
| What it reads | Residuals after the design is removed: components whose share of variance exceeds that of 20 row-wise permutations (p <= 0.10) | The canonical-correlation plot: the k with the largest gap between the all-feature curve and the control-feature curve, bounded (2 below 15 samples in proteomics, 3 in transcriptomics) |
| Needs | No control set; no missing values (run on the features complete in every sample) | A control set narrower than all features |
| Defined for | Every RUV option | Options with a separate control set (supplied list, first-pass nulls, Legacy); undefined with all features as controls (the two curves coincide) |
| Measured (KIT-124) | Mean 1.00 where the truth is 1; 3.00 where it is 3; 1.94 where it is 2; 0 in most no-factor data sets (S1, S2, S4), 1.35 in transcriptomics | 2.9 where the truth is 1; 1.55 where it is 3 (transcriptomics power 0.16 against 0.76) |
| Arms it ran in | 7, 8, 9, 10, 15, 16 | 3, 4, 5, 6 |
| For | Counts the factors; no control set needed; seeded and reproducible; the rule behind the recommended design's figures | The package's established diagnostic; ties k to the chosen controls; no new dependency; the rule behind Legacy's and first-pass nulls' figures |
| Against | Borrowed from surrogate variable analysis, not from an RUV paper; the design is removed first, so a factor that follows the groups can be under-counted; 20 permutations and a 0.10 threshold are coarse at 3 per group; optional dependency (`sva`) | Over-estimates at one factor, under-estimates at three; measures how well the controls avoid the biology, not how many factors there are; undefined for the recommended design |

Not yet measured, and needed before both rules can be selectable with a figure beside each: supplied controls with the separation rule's k (arm 8's design), and first-pass nulls with the estimated count (arm 5's design).

Provenance for the k rules:

| Claim | Source | Checked how |
|---|---|---|
| The permutation estimate is Buja and Eyuboglu's (1992), with Leek's (2011, Biometrics) asymptotic alternative | `sva::num.sv` help page and source, sva 3.60.0 | Read in the simulation library, 2026-11-18 |
| k is hard to estimate for RUV; `getK` (RUV-4 only) "will not be the best choice in all cases" and is "often a poor choice of K for use with RUV2" | `ruv::getK` help page, citing Gagnon-Bartsch and Speed 2012 (Biostatistics 13:539) and Gagnon-Bartsch, Jacob, Speed 2013 (Berkeley technical report 820) | Help page read 2026-11-18; the two papers themselves not read |
| RUVSeq's worked example uses k = 1 without a rule | RUVSeq vignette (section 1 above) | Read |
| RUV-III: no clear rule for k | Molania et al. 2019 (section 1 above) | Read |
| The canonical-correlation plot | `ruv::ruv_cancorplot`; the separation rule on it is MultiScholaR's own | Code read; no paper found that proposes the rule |

Citation check of 2026-11-18 (a separate reading pass; quotes copied from the fetched texts):

| Claim | Citation | Quote and place | Read as |
|---|---|---|---|
| Permutation version of parallel analysis | Buja A, Eyuboglu N (1992). Remarks on Parallel Analysis. Multivariate Behavioral Research 27(4):509-540. doi:10.1207/s15327906mbr2704_2 | "The idea of applying permutation test ideas to PA leads to a quasi-inferential, non-parametric version of PA" (abstract) | Abstract only |
| SVA counts factors by permuting each row of the residual matrix, crediting Buja and Eyuboglu | Leek JT, Storey JD (2007). Capturing Heterogeneity in Gene Expression Studies by Surrogate Variable Analysis. PLoS Genetics 3(9):e161. doi:10.1371/journal.pgen.0030161 | "Form a matrix R* by permuting each row of R independently to remove any structure in the matrix." (Methods, step 4); the significance level is "user-chosen" there, while `num.sv` fixes 0.10 and B = 20 | Full text |
| The "be" method is documented in the manual page, not in the 2012 package note | `sva::num.sv` help; Leek et al. (2012) Bioinformatics 28(6):882-883, doi:10.1093/bioinformatics/bts034, does not appear to describe it | "based on a permutation procedure originally prooposed by Buja and Eyuboglu 1992" [sic] | Help page read; the note's abstract only |
| The asymptotic alternative | Leek JT (2011). Biometrics 67(2):344-352. doi:10.1111/j.1541-0420.2010.01455.x | | Author manuscript |
| RUV-2: no clear way to find k, and an objection to SVA's testing approach | Gagnon-Bartsch JA, Speed TP (2012). Biostatistics 13(3):539-552. doi:10.1093/biostatistics/kxr034 | "there is no clear way to determine k. We recommend pursuing several approaches and exercising judgment."; "Using hypothesis testing to find k does not account for this bias-variance trade-off"; "This is the approach taken in SVA. However, we feel there are problems with this approach." | Preprint (Berkeley TR 800, sections 3.1 and 5); the first two and the first sentence of the third confirmed in the published text's index; published pages not confirmed |
| RUV-4: an automatic estimator exists (`ruv::getK`), hedged; RUV-inv avoids k | Gagnon-Bartsch JA, Jacob L, Speed TP (2013). Technical Report 820, Statistics, UC Berkeley | "We are unaware of any good algorithmic way to estimate k, and we feel there is an important role for human judgment." (3.6.6, p. 27); "we do not advise relying solely on k-hat" (5.4.3, p. 67); "both SVA and LEAPP estimate k using the method of Buja and Eyuboglu (1992)" (6, p. 71) | Full PDF |
| RUVSeq: k by judgment | Risso D, Ngai J, Speed TP, Dudoit S (2014). Nature Biotechnology 32(9):896-902. doi:10.1038/nbt.2931 | "The choice of k should be guided by considerations that include sample size, extent of technical effects captured by the first k factors, and extent of differential expression" (Online Methods) | Author manuscript |
| RUV-III with PRPS: a range of K, judged by metrics and biology | Molania R et al. (2023). Nature Biotechnology 41(1):82-95. doi:10.1038/s41587-022-01440-w | "RUV-III is generally robust to overestimating K but not always." (Methods) | Full text |
| RUV-III: a large k proposed | Jiang, Gagnon-Bartsch, Speed, arXiv:2310.01918; `ruv::RUVIII` help (k = NULL uses the maximum) | "we propose to choose k = K, a sufficiently large upper bound, such as K = m - s" (p. 6) | Passages on k only |
| Replicate-based RUV: judgment | Jacob L, Gagnon-Bartsch JA, Speed TP (2016). Biostatistics 17(1):16-28. doi:10.1093/biostatistics/kxv026 | "No single rule can be therefore given as to hyperparameter choice and judgment is necessary each time adjustments are performed" (5.1) | Full text |
| RUV-4 and SVA in one framework; factor count by bi-cross-validation | Wang J, Zhao Q, Hastie T, Owen AB (2017). Annals of Statistics 45(5):1863-1894. doi:10.1214/16-AOS1511 | "SVA [39] and RUV-4 [26] can be unified in this two-step statistical procedure." (section 1) | arXiv 1508.04178 |
| `num.sv` used to set k for RUVSeq in applied papers (precedent, not prescription) | Gatta et al. (2021) Schizophr Bull Open 2:sgab002, doi:10.1093/schizbullopen/sgab002; Bingen et al. (2022) Front Genet 13:1058817, doi:10.3389/fgene.2022.1058817; Marquez et al. (2025) FEMS Microbiol Ecol 101:fiaf053, doi:10.1093/femsec/fiaf053; Munzir et al. (2026) J Clin Med 15:4968, doi:10.3390/jcm15134968 | "The method of Buja and Eyuboglu implemented in the sva R package ... indicated that six latent variables could be used" (first of the four) | Full texts; two overlapping author groups; none gives a methodological reason |
| The separation rule on the canonical-correlation plot has no published source | `ruv_cancorplot` help: "Plots, as a function of k, the square of of the first canonical correlation of X and the first k left singular vectors of Y (and also, similarly, Y[,ctl])." [sic] | No rule in the `ruv` package, TR 820, arXiv 2310.01918 or a full-text search | Searched |

Reading for the paper: the estimated count is the SVA convention applied to RUV, with applied precedent and against the RUV-2 authors' stated reservation about test-based counts; its support here is the simulation, not a method paper. The RUV papers ask for judgment with diagnostics, which the user-chosen k and the diagnostic plots provide. The separation rule is MultiScholaR's own heuristic and is described as such. Quotes from preprints or author manuscripts (Gagnon-Bartsch and Speed 2012, Wang 2017, Risso 2014) are to be checked against the typeset articles before they are used. The CRAN `ruv` DESCRIPTION gives the 2019 NAR DOI for the 2012 paper: do not copy DOIs from it.

## 6. Decision record for the Methods section (the user's decisions, 2026-11-18)

One row per decision that the paper's Methods must state. "Basis" names the evidence or the convention; the measured figures are in `scored/` and the supplementary figures of KIT-144. Kept here because the tickets are not in the repository.

| # | Decision | Basis | State in the Methods as |
|---|---|---|---|
| M1 | Every RUV design is offered in every omic, in tiers (recommended; valid with conditions; legacy, warned; none), each with its measured false discovery proportion and note shown beside it | User: not the tool's place to decide for every user; section 4 | The option set and the tier of each option |
| M2 | Default: groups as pseudo-replicate sets, all features as negative controls, unwanted factors in the test model; pre-selected in all four omics | Sections 1 and 3: RUVSeq's RUVs usage; factors in the model is the RUVSeq and limma convention; arm 7 measured at the ceiling in three omics and within the limit in the fourth | The default, with the simulation figures; qualifiers from the per-data-set distributions (lipid/metab spread; classic proteomics at 12 per group) |
| M3 | Today's design is kept unchanged as "Legacy" in every omic, with its measured figure (0.30 to 0.56) | Reproducibility of earlier analyses; the departures of section 2 | A legacy option that reproduces earlier releases, and why it is not the default |
| M4 | DIA with limpa: the default replaces the limpa workflow script's control rule (least significant 1 to 30% by sweep); the sweep stays as Legacy. Figures measured on the `fastlimpadpc` engine, not on `limpa` | Section 3 (0.435 against 0.055); user: accepted provided the legacy sweep stays | The default for DIA; the engine the figures were measured on |
| M5 | Number of unwanted factors: the estimated count (`sva::num.sv`, method "be", B = 20, seeded) wherever defined; the package's separation rule selectable where a separate control set exists; a user-set number always possible; a diagnostic plot for either rule | Section 5: the count returned the true number in simulation; it is the SVA convention applied to RUV, with applied precedent and against the RUV-2 authors' reservation; the separation rule is MultiScholaR's own heuristic | Both rules, their provenance stated as here, and which rule each figure was measured with |
| M6 | Technical replicates: pre-selected as the replicate structure when a design declares them; transcriptomics sums technical replicates before analysis (edgeR User's Guide 4.9.1), so the option is shown disabled there with the reason; designs with one unit per group are refused in every omic | RUV-III's original use of technical replicates; standard RNA-seq practice; the void technical-duplicate cells before the repair | Handling of technical replicates per omic; the refusal |
| M7 | First-pass null controls offered in every omic, including DIA and lipid/metab where it measured poorly (0.149, 0.239), with that note | M1; RUVg's in-silico empirical controls | An option with stated conditions of use |
| M8 | Lipidomics and metabolomics: the level shift after RUV-III correction is repaired in every option including Legacy, so Legacy does not reproduce earlier releases bit for bit there (a constant per feature in the corrected table; contrasts unchanged on complete data) | A repaired defect is not an option; parity with proteomics, where it was repaired earlier | A correction noted in NEWS; a sentence if earlier results are compared |
| M9 | A feature unobserved in one whole group is still tested in the comparisons between groups where it is observed (not dropped from every comparison), flagged as not observed in that group; guarded by a test against the two-group fit | limma's default handling of missing values; today's behaviour; parity with proteomics; such features are often below detection, which is biology. Departs from the simulators' kernel, which dropped them (the simulations were mostly two-group) | Handling of features missing in a whole group |
| M10 | Under factors-in-the-model designs the test reads the uncorrected (normalised) data; displayed and exported intensities, heatmaps, PCA and the integration input are the RUV-corrected matrix. Fold changes and p-values come from the model, so a figure from the corrected matrix need not reproduce them exactly (difference up to 0.18 log2 units in the audit's data) | RUVSeq vignette (normalised counts for exploration, test on original data with W); limma `removeBatchEffect` (for plotting, not before modelling); one corrected matrix in the tool; parity with transcriptomics and DIA. The model-consistent alternative (data minus the fitted W terms) was considered and not taken | What is tested and what is displayed, stated explicitly; column headers and report text say so |
| M11 | A warning when an estimated unwanted factor follows the groups: threshold set on one set of simulated data and accepted on another, by a rule fixed beforehand (fires in at most 10% of unrelated runs and at least 80% of related ones); otherwise the value is shown with no verdict | KIT-129 section 3's declaration; arms fail mostly under confounding | The warning and how its threshold was set, or that none could be set |
| M12 | `fastduckdb` is a required dependency of the package and of every installer | User: the architecture depends on it; without it the lipid/metab intensity filter is refused | Software requirements; availability of the dependency |
| M13 | Simulations after 2026-11-18 run on GitHub-hosted runners; their results sit beside the laptop's only after an equivalence check under the fixed tolerance rule (relative 1e-6, absolute floor 1e-9, discrete values exact); the per-simulator time ceilings of the laptop runs do not apply to cells completed there (one attempt per cell) | User's platform decision; the rule already fixed for CI | Where the simulations ran and how equivalence was shown |
| M14 | Classic proteomics, factors-in-the-model designs: the test reads every protein of the normalised matrix (the population the simulations tested), not only the proteins that survive RUV-III-C and the missing-value clean-up of the corrected matrix; a protein tested without a corrected value is flagged, and has no displayed intensity | The measured figures are for the unfiltered population (audit of `kit-126`, F-03: 786 against 776 proteins in the centre scenario, 710 against 625 at 3 per group with k = 3); a filter computed on values the test does not read has no statistical basis for the test; consistent with M9. Decided by best practice at the user's request ("best practices here and most rigorous publishable outcome") | Which proteins are tested; that the tested set and the displayed set can differ by the proteins RUV-III-C cannot correct |

Still to be decided or measured before the Methods can state them: the three further arms (each selectable combination of control set and k rule), the technical-duplicate re-measurement, the warning's threshold, the extra-hard scenarios of three simulators, the real-data null check, the audits' open findings on the kernels.

## Sources

- RUVSeq vignette (Risso): https://rdrr.io/github/drisso/RUVSeq/f/vignettes/RUVSeq.Rmd
- limma `removeBatchEffect`: https://rdrr.io/bioc/limma/man/removeBatchEffect.html
- Nygaard, Rødland, Hovig 2016: https://pmc.ncbi.nlm.nih.gov/articles/PMC4679072
- Molania et al. 2019 (RUV-III): https://pmc.ncbi.nlm.nih.gov/articles/PMC6614807
- Jiang, Gagnon-Bartsch, Speed (theory of RUV-III): https://arxiv.org/abs/2310.01918
- Poulos et al. 2020 (RUV-III-C): https://pmc.ncbi.nlm.nih.gov/articles/PMC7393074
- Salim et al. 2022 (RUV-III-NB): https://pmc.ncbi.nlm.nih.gov/articles/PMC9458465
- Kim et al. 2021 (hRUV): https://pmc.ncbi.nlm.nih.gov/articles/PMC8371158

Full texts read through Europe PMC's REST service (PMC blocked automated fetches); quotes are short and verbatim. Not read: Molania et al. 2023 (PRPS, Nat Biotechnol), Gagnon-Bartsch and Speed 2012, Jacob et al. 2016, De Livera et al. 2015.
