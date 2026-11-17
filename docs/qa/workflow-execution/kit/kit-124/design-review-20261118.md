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
| Today's controls, W in the model (arm 4) | 0.322; 0.90 | 0.435; 0.90 | (arm 3) | 0.162; 0.17 |
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
| To verify before citing: Leek and Storey 2007 (PLoS Genet) and Leek et al. 2012 (Bioinformatics) for sva; Buja and Eyuboglu 1992, "Remarks on parallel analysis" (Multivariate Behavioral Research 27:509); Wang, Zhao, Hastie, Owen 2017 (Ann. Statist.) for RUV-4 and SVA in one framework with an estimated factor count | Recalled, not read this session | Not checked |

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
