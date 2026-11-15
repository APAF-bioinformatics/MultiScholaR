# KIT-124 W0: RUV convention, inventory and consumers

Reading only. No R code or test was changed. Code references are to branch `kit-123` (HEAD `5fe9e6dd`), which is
`performance-v2` plus two commits: `3e43d7a9` (R1, k bounded in the classic sweep) and `5fe9e6dd` (R2, levels restored
after correction). Where the two branches behave differently the cell says which is which.

Short forms: "controls" = negative controls; "W" = the estimated unwanted factors (samples by k); "level" = a feature's
mean over samples; "gap" = canonical correlation with the group of all features minus that of the controls, at a given
K, as `ruv::ruv_cancorplot` computes it.

## 0. How the sources were read

Each row says what was actually opened. Quotes are short and verbatim; where a fetch tool returned text in pieces the
quote was rejoined from consecutive pieces and should be spot-checked before it is printed anywhere else. Nothing in
this document is filled from memory.

| ID | Source | Read | Where |
|---|---|---|---|
| S1 | Poulos et al. 2020, Nat Commun 11:3793 (RUV-III-C) | Full text incl. Methods. Supplementary Notes not read. | https://europepmc.org/article/PMC/PMC7393074 |
| S2 | `RUVIIIC` 1.0.19 | R source, roxygen, tests, README. No vignette exists. Installed copy called. | https://github.com/cran/RUVIIIC ; https://github.com/CMRI-ProCan/RUV-III-C |
| S3 | Molania et al. 2019, NAR 47:6073 (RUV-III) | Full text. The "Supplementary Method" with the algebra not read. | https://europepmc.org/article/PMC/PMC6614807 |
| S4 | Molania et al. 2023, Nat Biotechnol 41:82 (RUV-III-PRPS) | Full text incl. Methods. Supplementary not read. | https://europepmc.org/article/PMC/PMC9849124 |
| S5 | The PRPS authors' code, `RUV_III_PRPS` | Read in full. | https://raw.githubusercontent.com/RMolania/TCGA_PanCancer_UnwantedVariation/master/Scripts/Libraries_HelperFunctions_ForTcgaReadVignette.R (lines 325-388) |
| S6 | Third-party `runRUV_III_PRPS` | Read in full. | https://raw.githubusercontent.com/AbhishekSinha28/tgcapkg/master/R/RUV_III.R |
| S7 | `ruv` package: `RUVIII`, `RUV1`, `ruv_cancorplot`, `getK`, `RUV4`, `RUVinv` | Source and help of 0.9.7.2 on the CRAN mirror; the installed 0.9.7.1 printed (`RUVIII`, `RUV1`, `ruv_cancorplot`). `RUV4`/`RUVinv` help read in part. | https://github.com/cran/ruv |
| S8 | Gagnon-Bartsch & Speed 2012, Biostatistics 13:539 (RUV-2) | Main text. Supplement not read. | https://pmc.ncbi.nlm.nih.gov/articles/PMC3577104/ |
| S9 | Gagnon-Bartsch, Jacob & Speed 2013, Berkeley Tech Report 820 (RUV-4, RUV-inv) | Full text. | https://statistics.berkeley.edu/sites/default/files/tech-reports/ruv.pdf |
| S10 | Risso et al. 2014, Nat Biotechnol 32:896 (RUVSeq) | Full text incl. Methods. | https://pmc.ncbi.nlm.nih.gov/articles/PMC4404308/ |
| S11 | RUVSeq vignette (1.46.0) and installed `RUVs` source | Read; `RUVs` printed. | https://bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.html |
| S12 | limma `removeBatchEffect` help and installed source (3.68.4) | Read. | https://rdrr.io/bioc/limma/man/removeBatchEffect.html |
| S13 | limpa vignette and reference manual (1.4.2), `dpcDE` | Read. | https://bioconductor.org/packages/release/bioc/vignettes/limpa/inst/doc/limpa.html ; https://bioconductor.org/packages/release/bioc/manuals/limpa/man/limpa.pdf |
| S14 | Nygaard, Rødland, Hovig 2016, Biostatistics 17:29 | Full text. | https://pmc.ncbi.nlm.nih.gov/articles/PMC4679072/ |
| S15 | Li, Zhang, Patil, Johnson 2023, Biostatistics 24:635 | Main text at the publisher, in pieces; one step of ComBat+Cor seen only as a paraphrase. | https://academic.oup.com/biostatistics/article/24/3/635/6514923 |
| S16 | Morris, White, Crowther 2019, Stat Med 38:2074 | Full text. Table 6's formulas came through garbled. | https://pmc.ncbi.nlm.nih.gov/articles/PMC6492164/ |
| S17 | Weber et al. 2019, Genome Biol 20:125 | Full text. | https://pmc.ncbi.nlm.nih.gov/articles/PMC6584985/ |
| S18 | Bradley 1978, BJMSP 31:144 | **Not read** (publisher refused). Secondary sources only, listed in section 1.8. | |
| S19 | Dunn et al. 2011, Nat Protoc 6:1060 | Full text (author-hosted PDF). | https://anivet.au.dk/fileadmin/user_upload/Dunn_Nature_Protocols_2011_6_1060_1083__untargeted_.pdf |
| S20 | Broadhurst et al. 2018, Metabolomics 14:72 | Full text. | https://pmc.ncbi.nlm.nih.gov/articles/PMC5960010/ |
| S21 | Fan et al. 2019, Anal Chem 91:3590 (SERRF) | Full text (author manuscript). | https://escholarship.org/content/qt59s2p8hw/qt59s2p8hw.pdf |
| S22 | Sysi-Aho et al. 2007, BMC Bioinformatics 8:93 (NOMIS) | Full text. | https://pmc.ncbi.nlm.nih.gov/articles/PMC1838434/ |
| S23 | De Livera et al. 2015, Anal Chem 87:3606 | Full text. | https://pmc.ncbi.nlm.nih.gov/articles/PMC4544854/ |
| S24 | Kim et al. 2021, Nat Commun 12:4992 (hRUV) | Full text. | https://pmc.ncbi.nlm.nih.gov/articles/PMC8371158/ |
| S25 | Köfeler et al. 2021, J Lipid Res 62:100138; Lipidomics Standards Initiative pages | Read. | https://pmc.ncbi.nlm.nih.gov/articles/PMC8585648/ ; https://lipidomicstandards.org/lipid-species-quantification/ |
| S26 | Märtens et al. 2023, Metabolites 13:665; Wehrens et al. 2016, Metabolomics 12:88 | Full text. | https://edoc.mdc-berlin.de/id/eprint/23395/1/23395oa.pdf ; https://pmc.ncbi.nlm.nih.gov/articles/PMC4796354/ |
| — | De Livera 2012 (Anal Chem 84:10768), Redestig 2009 (CCMN), Han & Li 2022 (Mass Spectrom Rev) | **Not read** (no open copy). | |

Installed versions printed from the pinned library: R 4.6.1, ruv 0.9.7.1, RUVIIIC 1.0.19, RUVSeq 1.46.0, limma 3.68.4,
limpa 1.4.2, edgeR 4.10.1.

## 1. Convention

### 1.1 What the references say, per choice

The sources do not agree with each other on three of the six choices. The disagreements are kept here and carried into
the verdicts.

**(a) Controls**

| Source | What it says or does |
|---|---|
| S1 (RUV-III-C paper) | Main dataset: "As negative controls, we selected peptides from yeast proteins (n = 1622) and scaled these such that their mean intensities were identical across samples." These are a constructed set from a mixed-in proteome, not synthetic spike-ins and not chosen by a statistic. Multi-lab dataset: of 8,555 human peptides detected in every sample, 500 were held out and "the remaining 8,055 were used as negative control peptides". General: "Endogenous genes will generally be a better choice of negative controls than will spike-ins"; "There are no pre-defined rules regarding the choice of negative controls, and domain-specific knowledge and human judgement are most important". |
| S2 (`RUVIIIC`) | `controls`: "The names of the control variables which are known to be constant across the observations". `RUVIII_C_Varying`: "Out of the potential negative control peptides, it selects those that are always non-missing across the selected subset." README example: all non-spiked human peptides found in every sample, 500 of them. |
| S3 (RUV-III) | "our approach to negative controls is pragmatic: if using a given gene in RUV-III as one of the set of negative control genes helps, as indicated by various measures, then whether or not it is an ideal negative control gene is moot"; "We will frequently get very good results using the entire set of genes being studied as negative controls". Used: the 500 or 340 genes of lowest variance, 9 of 236 of lowest variance, or all genes. |
| S4 (PRPS) | From the literature (housekeeping) or the data. One dataset: ANOVA against the biological subtype, "the genes with lowest F-statistics were selected (~1,000 genes)"; others add criteria that the gene shows the unwanted variation. |
| S8 (RUV-2) | "Negative control genes are genes known a priori not to be differentially expressed". On finding them in the data: "a researcher may wish to find genes whose expression levels are not highly correlated with the factor of interest, label these genes as negative controls, and then adjust via RUV-2. However, we feel this approach is misguided." Reason: "Excluding genes correlated with the factor of interest would bias our estimate of the unwanted factors." |
| S9 (RUV-4) | "Discovering empirical controls is often feasible, but may require some care." Example: "rank the genes by p-value and designate all but the top 100 as empirical negative controls". In its simulations: "all genes whose RUV-rinv FDR-adjusted p-values is greater that 0.5". "the strategy of empirical controls is much "safer" with RUV-4 than it is with RUV-2"; "We are unaware of any argument that would guarantee that the strategy of empirical controls will always lead to better results"; "it is important that the method used for the initial analysis exhibit good control of the type 1 error rate." |
| S10, S11 (RUVSeq) | "in silico" empirical controls: "the Jc least significantly DE genes based on a first-pass DE analysis". Used (an example): "all but the top 5,000 DE genes, as ranked by edgeR p-values" (15,839 and 16,500 genes in the paper; about three quarters of the filtered genes in the vignette, by arithmetic). For RUVs: "RUVs is robust to that choice and that using all the genes works well in practice". Restriction: empirical controls "are applicable only to classical DE settings". |
| S7 (`ruv`) | The help pages give no advice on choosing controls. |

No source selects the size of the control set by an objective, and none uses the gap for it. `ruv_cancorplot`'s help
says only: "Plots, as a function of k, the square of of the first canonical correlation of X and the first k left
singular vectors of Y (and also, similarly, Y[,ctl])." It gives no rule for reading controls or k off the plot. S9 does
not mention canonical correlation at all (searched).

**(b) Replicate structure**

| Source | What it says or does |
|---|---|
| S1 | Technical replicates: "we randomly partitioned the technical replicates of each of the eight sample types into groups of triplicates acquired on different instruments, and we then declared as technical replicates only those acquisitions corresponding to one of these groups." "RUV-III does not require that data from all samples be acquired with technical replicates." |
| S3 | Technical replicates and, in one dataset, "pseudo-replicate pairs". "If an nCounter gene expression dataset has very few or no replicates, RUV-III cannot be used." |
| S4 | Pseudo-replicates of pseudo-samples: averages of samples "roughly homogeneous with respect to the unwanted variation and biology", one per batch, grouped per biological subpopulation. The pseudo-samples are removed before downstream analysis. |
| S10, S11 (RUVs) | "replicate/negative control samples for which the covariates of interest are constant". Vignette example: "we can consider the three treated and the three control samples as replicate groups" (`makeGroups(x)`). Assumption: "both RUVr and RUVs assume that the unwanted factors are not correlated with the covariates of interest". |
| S24 (hRUV) | Replicates of the same sample: pooled QC, "short replicates" within a batch and "batch replicates" across batches. |

So "the group as the replicate set" is the RUVs vignette's example and is close to the PRPS idea, and it is not what
the RUV-III-C paper or the RUV-III paper did.

**(c) k**

| Source | What it says or does |
|---|---|
| S1 | "we chose the value of k to be 20, as this gave a high median Pearson correlation between all human peptide intensities and ovarian cancer tissue proportion across the dilution series"; multi-lab set k = 11. In general k is set "using information such as measurements on positive control peptides (when differential abundance is the goal), graphical displays or by exploiting some special context-dependent feature of the study." README: "Non-dilution datasets will likely use a smaller value of `k`." |
| S3 | "there are no clear rules for determining k"; repeat with several k and judge by metrics and prior knowledge. "RUV-III is generally robust to overestimating k, but not always." Used k = 1, 6, 10 and the maximum. |
| S4, S5 | Same approach; the authors' vignette code uses k = 20 and k = 12. |
| S8 | "there is no clear way to determine k"; RLE plots, p-value histograms, positive controls. "Using hypothesis testing to find k does not account for this bias-variance trade-off". |
| S9 | "We are unaware of any good algorithmic way to estimate k, and we feel there is an important role for human judgment." RUV-inv "eliminates any need to estimate the number of unwanted factors". `getK` is documented for RUV-4 only, with a warning that it "will not be the best choice in all cases". |
| S10, S11 | k "should be guided by considerations that include sample size"; "the small sample size (n = 6) for the Zebrafish dataset only allows the estimation of one or two factors"; every vignette example uses k = 1. |
| S7 | `RUVIII`'s default is `k = NULL`, "the maximum possible value of k"; `k = 0` is allowed and makes no adjustment. |
| S24 | k fixed: "defaulting to 5". |

No source gives an algorithm, and none bounds k by a table of sample sizes. The published values are not small in the
RUV-III family (20, 11, 12, 10, 6, 5) and are 1 in RUVSeq's examples.

**(d) The test, and (e) what the corrected matrix is for**

| Source | What it says or does |
|---|---|
| S3 (RUV-III) | Tests on corrected values: "Differential expression analyses of all normalized Nanostring datasets in this paper were performed using the limma R package. The design matrix was created using the model.matrix function for only biological factor of interests." |
| S1 (RUV-III-C) | Corrected values (after the paper's own replacement of missing values) are used for coefficients of variation, correlation with the dilution series, unpaired t-tests after averaging each set of technical triplicates, and machine learning. No limma-type model on corrected values, and no W in a model. |
| S4 (PRPS) | Corrected values used for subtype separation, co-expression and survival. Its tests on corrected values (Wilcoxon) compare batches, as a diagnostic. |
| S8, S9 (RUV-2, RUV-4) | W in the model: "perform factor analysis on just the negative control genes and incorporate the resulting factors into a linear regression model"; RUV-4 "includes the estimated unwanted factors as covariates in a regression model". |
| S10, S11 (RUVSeq) | W in the model, original counts: "RUV normalization is performed by including the estimated W matrix in the GLM." Vignette, of the normalised counts: "These counts should be used only for exploration." and "It is important that subsequent DE analysis be done on the original counts". Paper: "removing Wα from the original counts bears the risk of removing part of the effect of X". |
| S12 (limma) | "This function is not intended to be used prior to linear modelling. For linear modelling, it is better to include the batch factors in the linear model." |
| S14 (Nygaard) | Adjusting for batch while keeping group differences "may systematically induce incorrect group differences in downstream analyses when groups are distributed between the batches in an unbalanced manner". Balanced: "results should be reliable". Advice: "account for batch in the statistical analysis." |
| S15 (Li) | "two-step approaches introduce a correlation structure in the corrected data, which, if ignored, may lead to either exaggerated or diminished significance". "the endogeneity issue does not exist in balanced batch-group design." It still recommends two-step ComBat for balanced designs and proposes ComBat+Cor (generalised least squares on the corrected data) for unbalanced ones. |
| S13 (limpa) | "limpa works with any design matrix, with any combination of explanatory factors and covariates." "limpa differential expression analyses must be done using `dpcDE`". Neither the vignette nor the manual mentions RUV, batch or unwanted variation (searched). |

The disagreement: the RUV-III papers test on corrected values with a design of the biological factors only (S3 does
exactly that with limma); the RUV-2, RUV-4 and RUVSeq papers, limma's help and S14 put the factors in the model and
keep corrected values for display. None of S1, S3, S4 discusses the degrees of freedom or the correlation that S14 and
S15 describe.

**(f) Level of each feature**

| Source | What its code does |
|---|---|
| S7 `ruv::RUVIII` | No centring: `W = Y[, ctl] %*% t(ac) %*% solve(ac %*% t(ac))`, `newY = Y - W %*% alpha`. With the default `eta = NULL`, `RUV1` returns Y unchanged. |
| S2 `RUVIIIC` | No centring in the R prototype or the C++ ("No additional transformation is applied to the input data matrix"). Its own test asserts equality with default `ruv::RUVIII`. |
| S5 PRPS authors' code | Centres the controls for W only: `mu <- colMeans(Y)`, `Y_stand <- Y - mu_mat`, `W <- Y_stand[, ctl] %*% t(ac) %*% solve(ac %*% t(ac))`, `newY <- Y - W %*% alpha`; `alpha` from uncentred `Y`. The paper states it: model "Y=1μ+MXβ+Wα+ε" and "Estimate W by regressing the centered negative controls". |
| S6 third-party copy | Byte-identical to `tcgaCleaneR/R/RUV_III.R` in the authors' repository; same centring lines as S5. **The authors' own code does centre the controls; the third-party copy did not add it.** |
| S1, S3 | Say nothing about centring or levels. S3's algebra is in a supplement that was not read. |
| S12 limma | `removeBatchEffect` centres its covariates (`covariates <- t(t(covariates) - colMeans(covariates))`), so levels are kept. |
| S11 `RUVSeq::RUVs` | `W <- Y[, cIdx] %*% ...` on uncentred log counts; `correctedY <- Y - W %*% a`. Not centred. |

Checked here on a toy matrix (12 samples, 200 features at about 20, 40 controls, k = 2, two groups as replicates):
`ruv::RUVIII` and `RUVIIIC::RUVIII_C_Varying` agree to 2e-11 and both move the feature means by 2.1 on average. So the
two reference packages shift levels; only the PRPS code keeps them. Derived, not quoted: with centred controls every
column of W sums to zero, so each feature's mean is unchanged for a positive k.

### 1.2 Proteomics classic (DIA with IQ, LFQ, TMT)

Run by `diaNormalizationOperation` (`R/utils_workflow_dia_products.R:337-380`): `resolveProtRuvParameters`
(`R/utils_workflow_dia_norm.R:53-102`), then `ruvIII_C_Varying` (`R/func_prot_s4_norm_methods.R:204-288`).

| Choice | What the code does | Reference | Verdict |
|---|---|---|---|
| (a) Controls | Per-protein one-way ANOVA on the group (`func_prot_norm_optimization_helpers.R:557-563`), groups whose name contains "pool" left out (:500); q-values, `qvalue` in the automatic sweep (:690), BH in manual mode (`utils_workflow_dia_norm.R:89`); proteins with q > 0.05 kept and the N of highest q taken (:582-588); N = the percentage of all proteins. The sweep tries 1 to 20% (:686; the GUI allows up to 50, `utils_omic_gui_config.R:1234`) and keeps the percentage with the highest composite score: the largest gap over K, times one minus a penalty that grows with k (:172-199, :219-244, :923). | S1: constructed or all-complete controls, "human judgement". S9, S10: empirical controls are wide sets (all but the top 100; q > 0.5; all but the top 5,000). S8: choosing features uncorrelated with the factor of interest is "misguided". No source chooses the size by an objective. | **Differ.** A narrow set of the least significant, sized by an objective no source defines and which rewards controls that are uncorrelated with the group. |
| (b) Replicates | Technical replicates when the design declares them and they are not spanned by the group, else the RUV grouping variable (`resolveProtRuvReplicateBlocks`, :417-433; rule in `resolveProtTechnicalReplicateBlock`, :399-415). The structure used is recorded (`func_prot_s4_norm_methods.R:245-248`). | S1, S3: technical replicates. S11: the group, as an example for RUVs. S4: pseudo-replicates within a biological group. | Technical replicates: **agree** with S1. Group fallback: **differ** from S1 and S3, same as the S11 example. |
| (c) k | Elbow: the smallest K whose gap is within 5% of the largest gap; 1 when the largest gap is under 0.05 (`findBestKElbow`, :64-91). **kit-123:** the curve is cut at `max_acceptable_k` (2, 3, 4, 5 for under 15, 40, 80 and more samples, :270-289) before the score and the elbow (:885-894), so k is 1 to 5. **performance-v2:** the elbow read the whole curve; the limit entered only the penalty. Both: a percentage whose k exceeds its controls observed in every sample is skipped (:910). Manual mode: any k of 1 or more, default 3 (`utils_workflow_dia_norm.R:82-83`). | No rule in any source; judgement with positive controls and diagnostics (S1, S3, S8, S9). S1 used 20 and 11. | **Differ** (an algorithm where the sources ask for judgement). Whether it does better or worse is **unknown** until W3. |
| (d) Test | limma on the corrected table with the user's formula, default `~ 0 + group`; no W; `duplicateCorrelation` when technical replicates exist (`func_prot_da_model_methods.R:408-417`; `func_prot_da_model_stats.R:231-266`). Residual degrees of freedom are not reduced by k. | S3: the same (limma on corrected values, biological factors only). S8-S12, S14: factors in the model. | **Agree** with the RUV-III paper's practice; **differ** from RUV-2, RUV-4, RUVSeq, limma's help and S14. |
| (e) Corrected matrix | It replaces the protein table (`func_prot_s4_norm_methods.R:282`) and is the input of everything after it: section 4. | S1, S3, S4: corrected values feed downstream analysis. S11: exploration only. | **Agree** with S1/S3/S4; **differ** from S11. |
| (f) Level | **kit-123:** each protein's mean over the samples it was corrected in is set back to its mean before (`restoreRuvFeatureLevels`, :1152-1159, called at `func_prot_s4_norm_methods.R:269`). **performance-v2:** as `RUVIII_C_Varying` returns it, shifted. | S2, S7: shifted. S5: kept. | kit-123 **differs** from `RUVIIIC` by a constant per protein and **agrees** with S5 on complete data (by algebra; not tested here). performance-v2 **agrees** with `RUVIIIC`. |
| Level of application | Protein table. | S1 applies RUV-III-C to peptides ("we normalise the data on any given peptide p"). | **Differ.** |

### 1.3 DIA_LIMPA

Run by `diaRuvCovariatesOperation` (`R/utils_workflow_dia_products.R:782-853`) after the correlation filter:
`findBestNegCtrlPercentageCapped` (`func_prot_norm_optimization_helpers.R:1060-1142`), then `computeRuvCovariates`
(:1174-1205). Parameters: `resolveDIARuvCovariatesParameters` (`utils_workflow_dia_products.R:859-877`).

| Choice | What the code does | Reference | Verdict |
|---|---|---|---|
| (a) Controls | The same ANOVA selection, BH (:1065), 1 to 30% (:1061), on the complete DPC-quantified matrix. Score: the largest gap over K up to `min(max_acceptable_k, controls - 1)` (:1105-1107), with the same penalty; the first maximum wins (:1135). | As 1.2. limpa (S13) says nothing on RUV; the rule comes from the limpa script that KIT-041 brought in. | **Differ**, as 1.2. |
| (b) Replicates | The same rule as classic (:1187-1191). | As 1.2. | As 1.2. |
| (c) k | The K at the largest gap on the capped curve (:1110-1113): 1 to the cap (2 to 5 by sample size), and below the number of controls. | No rule in any source. | **Differ**; effect **unknown** until W3. |
| (d) Test | Covariates mode (default): `limpa::dpcDE` on the uncorrected EList with the user's formula plus `W1..Wk` (`func_prot_da_model_methods.R:378-392`). Correction mode: `dpcDE` on the EList whose `E` was replaced by `Y - W alpha`, plain formula, standard errors and observation counts kept as quantified (`utils_workflow_dia_products.R:817-833`); a design that still carries W columns is refused (`func_prot_da_model_methods.R:380-381`). | Covariates: the form of S8-S11, with W from RUV-III (as `RUVs` computes it) instead of RUV-4. S13 allows any covariates. Correction: the form of S3; S13 does not describe changing `E` while keeping the DPC standard errors. | Covariates: **agree** in form with S10/S11. Correction: **agree** with S3; its effect on `dpcDE`'s weights is **unknown**. |
| (e) Corrected matrix | Covariates mode: a QC matrix only (`ruv-adjusted-science`, :815; file `ruv_adjusted_protein_matrix.tsv`); the state keeps the quantified matrix and gains W in the design (:801). Correction mode: the state's table and the EList. | S11: exploration only. | Covariates: **agree** with S11. Correction: as 1.2 (e). |
| (f) Level | **kit-123:** W from the controls centred on their means (:1198-1199), `adjusted = Y - W alpha` (:1202): levels kept, in the QC matrix and in correction mode. **performance-v2:** uncentred W and `ruv::RUVIII`'s `newY`: shifted. In covariates mode the test is the same on both branches as long as the formula has an intercept or `0 + group` (derived: the constant part of an uncentred W lies in the span of the group columns). | S7: shifted. S5: kept. | kit-123 **differs** from `ruv::RUVIII`'s `newY` and W, **agrees** with S5. |

### 1.4 Transcriptomics

`trxRuvGrid` and `trxRuvs` (`R/func_trx_ruv.R:16-157`), `trxRuviii` (:165-191); operations in
`R/utils_workflow_trx_declarations.R` (:217-282, :394-405); the test at :290-340.

| Choice | What the code does | Reference | Verdict |
|---|---|---|---|
| (a) Controls | A first-pass edgeR quasi-likelihood F-test over all group coefficients (`func_trx_ruv.R:43-49`); the least significant percentage, at least 10 genes, no q cut-off (:60-62); percentages 5 to 30 by 5 (:18). Score: the largest gap over K up to the cap, times one minus a penalty of weight 0.05 (:66-76); first maximum wins (:119). | S10, S11: first-pass edgeR test, "all but the top 5,000" (an example; about three quarters of the genes). | Same idea; **differ** in size (5 to 30% against about 75%) and in sizing the set by an objective. |
| (b) Replicates | `RUVSeq::makeGroups(design$group)` (:146): the group. | S11's example. | **Agree.** |
| (c) k | The K at the largest gap up to the cap: 3 below 15 samples, else 4 (:38, :82-90). Manual k above the cap is refused (`utils_workflow_trx_declarations.R:245-248`). | S10: one or two factors for n = 6; k = 1 in every example. | **Differ** (chosen from the data; up to 3 with fewer than 15 samples). |
| (d) Test | edgeR on the original counts with `~ 0 + group + W_1 + ... + W_k`, robust dispersion, quasi-likelihood F-test (:298-340). A design where W is confounded with group is refused (:315-317). | S11: W in the design, original counts (its example uses `glmLRT`). | **Agree** (quasi-likelihood F instead of the likelihood-ratio test of the example). |
| (e) Corrected matrix | Two. The display matrix: `limma::removeBatchEffect` of the TMM logCPM with W as covariates and `~ group` protected (`func_trx_ruv.R:153-154`), for figures and `ruv_normalised_log_values.tsv`. The RUVIII export: `ruv::RUVIII` on the logCPM with the same controls and k, M from the group or from technical replicates (:165-178), written for integration. RUVSeq's own normalised counts are returned by the kernel but not stored by the operation (`utils_workflow_trx_declarations.R:278`). | S11: normalised counts "only for exploration". S3/S4: RUV-III values feed downstream analysis. | Display: **agree** with S11. Export: a second correction that RUVSeq does not describe; **agrees** in use with S3/S4. |
| (f) Level | Display matrix: kept (limma centres the covariates). RUVIII export: `ruv::RUVIII`'s `newY`, not restored on either branch: shifted. | S12: kept. S7: shifted. | Display **agree**; export **agrees** with `ruv::RUVIII` and so shifts. |

### 1.5 Lipidomics and metabolomics

The two S4 files are the same text with the omic's name changed (`R/func_lipid_s4_ruv_methods.R`,
`R/func_metab_s4_ruv_methods.R`, compared here). Run by `omicsRuvOperation`
(`R/utils_workflow_omics_declarations.R:305-329`). Sweeps: `runLipidPerAssayRuvOptimization`
(`R/func_lipid_norm_ruv_optimization_helpers.R:33`), `runPerAssayRuvOptimization` (`R/func_metab_norm.R:353`). Each
assay (for example positive and negative mode) is treated separately.

| Choice | What the code does | Reference | Verdict |
|---|---|---|---|
| (a) Controls | The proteomics ANOVA helper per assay (`func_lipid_s4_ruv_control_methods.R:334`; `func_metab_s4_design_control_methods.R:417`), BH, q > 0.05, 1 to 20% (lipid helpers :49-50, :102; metab :369-370, :422). Score as classic, on the whole curve (lipid :250-296; metab :570-616). | Field: internal standards or quality-control metabolites (S23); all metabolites (S24); see section 2. S23: empirical controls "need to be found empirically" when standards are too few, with "Care needs to be taken that the empirical quality control metabolites are not influenced by the factors of interest." | **Differ** from S24; the S23 fallback in kind, sized by an objective S23 does not define. |
| (b) Replicates | The RUV grouping variable only (`func_lipid_s4_ruv_methods.R:690-696`; metab the same lines). The `technical_replicate_id` slot is read and not used (:413). | S24: replicates of the same sample. S1: technical replicates. | **Differ.** |
| (c) k | Elbow on the whole curve, no upper limit (lipid helpers :255; metab :575); `max_acceptable_k` enters only the penalty (:296; :616). Skipped when k exceeds the complete controls. Manual default 2 (:56; :376). KIT-123 R1 did not change these: both branches behave the same. | S24: fixed at 5. S23: "there is no "optimal" λ and k"; "chosen using context-specific knowledge". | **Differ**; effect **unknown** until W3. |
| (d) Test | limma `lmFit` on the corrected assay with the user's formula, default `~ 0 + group`; no W, no replicate block (`func_metab_da_core_engine.R:6-37`; `func_lipid_da_model.R:6-14`). | S23: RUV-2, a linear model with the factors in it, for differential abundance; RUV-random gives a normalised matrix for any use. S3: tests on corrected values. | **Agree** with S3 and RUV-random's use; **differ** from RUV-2. |
| (e) Corrected matrix | Replaces the assay data (:913) and feeds everything after it: section 4. | S24, S23 (RUV-random): a normalised matrix for downstream use. | **Agree.** |
| (f) Level | As `RUVIII_C_Varying` returns it; `restoreRuvFeatureLevels` is not called (its only caller is the proteomics method). Shifted on both branches. Expected from the code and the toy check in 1.1; not measured on lipid or metabolite data. | S2: shifted. S5: kept. | **Agree** with `RUVIIIC`; **differ** from S5 and from classic proteomics on kit-123. |

### 1.6 Peptide-level method

`ruvIII_C_Varying` for `PeptideQuantitativeData` (`R/func_pept_s4_norm_methods.R:841-918`) and
`findBestNegCtrlPercentagePeptides` (:235). No GUI or worker route calls them (section 3).

| Choice | What the code does | Reference | Verdict |
|---|---|---|---|
| (a) Controls | The ANOVA helper on the peptide matrix (:198); sweep 1 to 20% (:238), score as classic. | S1: constructed or all-complete peptides. | **Differ.** |
| (b) Replicates | The RUV grouping variable only (:880-882); no technical-replicate rule. | S1: technical replicates. | **Differ.** |
| (c) k | Elbow on the whole curve (:428); neither KIT-123 fix reaches it. | No rule. | **Differ.** |
| (d) Test | None of its own; the workbooks continue to protein quantification. | — | not applicable |
| (e) Corrected matrix | Replaces the peptide table and matrix (:913-914). | S1: the corrected peptide matrix is the product. | **Agree.** |
| (f) Level | As `RUVIII_C_Varying` returns it: shifted. | S2: shifted. | **Agree** with `RUVIIIC`. |

It is the only implementation at the level the RUV-III-C paper works at.

### 1.7 What all five share that no source has

- The control set is sized by an objective computed with the group: the gap. S8 warns against choosing controls for low
  correlation with the factor of interest; S9 and S10 use wide sets and fix their size in advance.
- k is read from the same curve by a rule (elbow, or largest gap) under a cap set by sample size.
- `ruv_cancorplot` takes its singular vectors from the uncentred matrix (`W1 = svd(Y)$u`), so on log-scale data the
  first vector mostly carries the overall level. No source discusses reading k from this plot.

### 1.8 The four standards the ticket names

| Standard | What the source says | The ticket's use |
|---|---|---|
| Morris, White, Crowther 2019 (S16) | "defining aims, data‐generating mechanisms, estimands, methods, and performance measures ("ADEMP")". The E is estimands; section 3.3 is "Estimands and other targets" and allows a null hypothesis as the target ("if the target is a null hypothesis, power and type I error rates will be of primary interest"). "we should present estimates of uncertainty (quantified as the Monte Carlo standard error". Repetitions: "key performance measures need to be estimated to an acceptable degree of precision"; worked example 95 × 5 / 0.5² = 1900. | **Correct.** "Targets" is the paper's wider word and fits a test. Setting the repetitions from a pilot's between-run standard deviation and a required standard error is what section 5.3 describes. The bounds of 200 and 1,000 repetitions are the ticket's own. |
| Weber et al. 2019 (S17) | "A neutral benchmark should include all available methods"; for a method's own evaluation, "the current best-performing methods (if known), a simple 'baseline' method". On measures: "The choice of metric depends on the type of method and data"; "subjectivity in the choice of datasets or evaluation metrics could bias the results." Simulated data give "a known true signal"; "Experimental data often do not contain a ground truth". | Current methods in the comparison: **correct**. "Measures are fixed in advance": **not in the text**. The words "in advance" occur once, about parameters; there is no "pre-specified" or "a priori". Weber supports only that the choice of measures can bias. Pre-specification is the ticket's own rule (and consistent with S16's planning), not Weber's. Simulated and real data: Weber describes both and their trade-off; a sentence requiring both was not found in what was read. |
| Bradley 1978 (S18) | Original **not read**. Secondary: Blanca et al. 2017, "According to Bradley's (1978) liberal criterion a statistical test is considered robust if the empirical Type I error rate is between .025 and .075 for a nominal alpha level of .05." (https://www.psicothema.com/pdf/4434.pdf). Zumbo & Jennings 2002 give the same band as the "very liberal criterion" and 0.045 to 0.055 as "fairly stringent" (https://www.redalyc.org/pdf/169/16923209.pdf). | The numbers match: 0.075 is 1.5 times 0.05. Two limits on the use. The criterion is two-sided and W3 rule 1 uses only the upper limit, so a test with a rate under 0.025 passes it. And it was proposed for a test's Type I error under violated assumptions; applying it to a whole pipeline is by analogy and should be said so. The name of the band varies between the secondary sources, so Bradley's wording should not be put in quotation marks. |
| Pre-specification | The ticket's own rule. | No source to check. |

## 2. Lipidomics and metabolomics convention

| Family | Needs | Does | Source |
|---|---|---|---|
| Pooled QC, LOESS (QC-RLSC) | A pooled QC sample injected through the run; injection order. | Per feature, a LOESS curve through the QC injections, then the whole run divided by it. "one QC sample at every fifth injection"; features with QC RSD over 20% removed. | S19 |
| Pooled QC, guidance | The same. | "it is recommended that a pooled QC is analysed at least every 10th sample." Corrections without pooled QC: "these methods are not recommended by the authors." | S20 |
| Pooled QC, SERRF | Pooled QC, injection order, batch. | A random forest per feature trained on the QC injections, with "the injection order, batch effect, and the intensity of the QCs of the other metabolites as predictors". "QC-based normalization approaches are becoming more popular." | S21 |
| Internal standards | Labelled or exogenous standards added before extraction. | Lipidomics: "at least one exogenous IS should be added per lipid class as a minimum requirement"; "Data should be normalized to appropriate internal standard and the sample amount". NOMIS: an optimal factor per feature from several standards. Untargeted metabolomics: "it is impossible to use internal standards comprehensively" (S20). | S25, S22, S20 |
| RUV with standards as controls | Internal standards or other "quality control metabolites". | RUV-2 for differential abundance only ("RUV-2 cannot be used prior to unsupervised analyses"); RUV-random gives "The normalized data matrix" for any use. | S23 (De Livera 2012 not read) |
| Replicate-based, hRUV | Pooled QC, replicates within a batch and replicates carried to the next batch; run order. | A smoother against run index within the batch, then RUV-III across batches on the replicates. "All metabolites were used as the negative controls for RUV" and k "was taken to be 5". | S24 |
| RUV with QC samples | Pooled QC. | Unwanted variation estimated from the QC samples; k = 3. "Probably because of the low number of QCs, also RUV is not able to arrive at the same quality level." | S26 (Wehrens) |

**Where RUV with data-chosen controls stands.** In what was read:

- The default for drift and batch in untargeted LC-MS is pooled-QC regression against injection order (S19, S20, S21,
  S26 Märtens). Lipidomics guidance asks for internal standards per class and is silent on statistical batch
  correction (S25).
- Data-chosen controls appear as a fallback inside the RUV literature: S23 allows them when standards are too few, by
  nearness to the standards or by "p-value and/or fold-change rankings", with the caution quoted in 1.5.
- hRUV tried controls chosen by p-value against all metabolites and "saw little difference"; it attributes the
  correction to the replicates. That is RUV-III with true sample replicates, not with the group as the replicate set.
- No guideline or independent comparison that was read recommends or ranks RUV with data-chosen controls and the group
  as replicates, which is what this package does for lipids and metabolites. The review most likely to (Han & Li 2022)
  could not be opened, so this is absence in what was read, not a negative finding.

**What the package offers beside RUV.**

- Internal standards: `normaliseUntransformedData(method = "ITSD")` for both omics
  (`R/func_lipid_s4_normalization_methods.R:1-69`; metabolomics at `R/func_metab_s4_transformation_methods.R:220-224`).
  Standards are found by a pattern (`(?i)itsd|istd|internal.?standard`, `R/func_lipid_norm_support_helpers.R:24`,
  `R/func_metab_norm_support_helpers.R:116`) or picked in a table. One factor per sample: the standards are combined
  (sum, mean or median) and each sample is scaled to the average. That is one global standard signal, not one standard
  per lipid class as S25 asks. The standards are removed after normalisation by default
  (`remove_itsd_after_norm = TRUE`), so they are not available to RUV as controls.
- Between-sample normalisation and a QC tab for the standards (`R/mod_lipid_qc_itsd.R`,
  `R/mod_metab_qc_itsd_server.R`).
- No pooled-QC drift correction (no LOESS on QC, no SERRF) and no use of injection order in normalisation: none found
  under `R/` (searched `loess`, `serrf`, `injection`, `run_order`; `run_order` in the design configuration is the
  order of samples in the design table, `R/func_omics_design_shared.R:466-478`).

**Can a design mark pooled QC samples?** Not as a role. The design table has `group`, `factor1` to `factor3`, `batch`,
`replicates` and `tech_reps` (`R/func_omics_design_shared.R:505-506`); there is no sample-type column. Pooled samples
are recognised by name only: a group whose name contains "pool" is left out of the ANOVA that picks controls
(`func_prot_norm_optimization_helpers.R:500`) and a group matching "pool" or "qc" is left out of the correlation
filter (`R/func_prot_s4_qc_correlation_methods.R:192`; lipid and metab use the pattern `pool`,
`R/func_lipid_s4_correlation_methods.R:265`, `R/func_metab_s4_qc_methods.R:742`). In RUV itself a pool group is one
more replicate set like any other group; nothing uses pooled QC injections to estimate drift.

## 3. Inventory

| Omic | Known batch in the test's model | Technical replicates and pooled QC in the design; how RUV uses them | Peptide-level RUV |
|---|---|---|---|
| Proteomics classic (DIA, LFQ, TMT) | Yes, by typing it. The formula is an editable text input, default `~ 0 + group` (`R/utils_omic_gui_config.R:1368-1369`); the design has a `batch` column (`func_omics_design_shared.R:505`); the fit builds the model from the string (`func_prot_da_model_stats.R:438`). Not exercised here. There is no batch selector, and the RUV tab's "skip" help points to this route ("Consider including batch as a covariate", :1244-1245). | `tech_reps` in the design sets `tech_rep_group` = group and replicate number (`func_omics_design_shared.R:651`, :687), which the object declares as `technical_replicate_id`. RUV: these blocks are M when at least one block has two runs and the blocks are not the groups themselves; else the group (`func_prot_norm_optimization_helpers.R:399-433`). The same block goes to `duplicateCorrelation` in the test (`func_prot_da_model_stats.R:231-252`). Pool groups: left out of control selection only. | Not called. |
| DIA_LIMPA | Yes, the same formula input; W terms are appended to it in covariates mode (`func_prot_da_model_methods.R:382`). | The same rule for M (`:1187-1191`) and for `dpcDE`'s block (`func_prot_da_model_methods.R:388-392`). | Not on the GUI route. The three DIA_LIMPA workbooks call it (below). |
| Transcriptomics | No. The model is shown, not edited (`utils_omic_gui_config.R:2129-2133`); the formula is always `~ 0 + group + W_1..W_k` (`utils_workflow_trx_declarations.R:298-301`) and any other is refused (:346-350). | `technical_replicate_id` exists on the class. RUVs always uses the group. Only the RUVIII export can use technical replicates for M, and only when a set has two or more runs within one group (`func_trx_ruv.R:180-191`; choice at `utils_omic_gui_config.R:2098-2101`). The edgeR test has no replicate block. | Not applicable. |
| Lipidomics | Yes, by typing it: editable formula, default `~ 0 + group` (`utils_omic_gui_config.R:1146-1147`); model built from the string (`func_lipid_da_model.R:14`). Not exercised here. | The slot exists (`func_lipid_s4_assay_data_class.R:57`); RUV ignores it and uses the group (`func_lipid_s4_ruv_methods.R:690-696`); the test has no block. Pool groups: left out of control selection. | Not applicable. |
| Metabolomics | As lipidomics (`func_metab_da_core_engine.R:13-15`). | As lipidomics (`func_metab_s4_objects.R:78`; `func_metab_s4_ruv_methods.R:690-696`). | Not applicable. |

Two further points from the GUI configuration. The RUV grouping variable is a choice on the tab, and its help says
"For TMT workflows with batches, select 'batch' to remove batch effects" (`utils_omic_gui_config.R:1228-1230`). With
`batch` as the grouping variable the controls are the proteins that differ least between batches and M treats the runs
of a batch as replicates, which protects batch differences instead of removing them; the same help text calls the
variable the "Experimental factor to preserve". No source describes this use. Second, no omic offers the known batch
as a model term through a control; where it is possible it is by editing the formula.

**Callers of the peptide-level method** (`ruvIII_C_Varying` for `PeptideQuantitativeData`,
`findBestNegCtrlPercentagePeptides`, `getNegCtrlProtAnovaPeptides`, `ruvCancorFast`):

- No operation declaration, module or worker step refers to them.
- Three workbooks call them: `inst/workbooks/proteomics/DIA_workflow_limpa.rmd:766, 783`,
  `DIA_workflow_limpa_starter.rmd:716, 733`, `DIA_workflow_limpa_experienced.rmd:761, 778`.
- Otherwise they are reachable only as exported functions. `ruvCancorFast` has no caller in `R/` or `inst/`.

Placement for the ticket: the peptide-level method is outside every GUI route. It can be scoped out of the simulations
without leaving a route untested; the three workbooks are the only users.

## 4. Consumers

Traced in `R/` and `inst/` on kit-123. "Corrected" means the RUV-corrected values; "uncorrected" the normalised values
before RUV.

### 4.1 Producers

| Route | Producer | What it writes |
|---|---|---|
| Proteomics classic | `ruvIII_C_Varying` (`func_prot_s4_norm_methods.R:204-288`) in `diaNormalizationOperation` (`utils_workflow_dia_products.R:337-380`) | Checkpoints `normalized-science`, `ruv-science`; state `ruv_corrected` (`normalized` when skipped, :595-609); `norm_data$ruv_normalized_obj` (:622). Files (:641-662): `normalized_protein_matrix_pre_ruv.tsv`, `ruv_optimization_results.RDS`, `ruv_corrected_{pca,density,rle,correlation,cancor}.png`, the composite figure; after the correlation filter `ruv_normalised_results_cln_with_replicates.tsv/.RDS`. |
| DIA_LIMPA | `diaRuvCovariatesOperation` (`utils_workflow_dia_products.R:782-853`) | `W1..Wk` in the design (:801); `@args$ruvCovariates` (:803-808); checkpoints `ruv-covariates-science`, `ruv-adjusted-science` (never a state), and in correction mode `ruv-corrected-science` with the table and the EList's `E` replaced (:817-833). State `ruv_covariates` or `ruv_corrected` (:933-953). `workflow$ruv_normalised_for_da_analysis_obj` (:972). Files (:977-987): `ruv_covariates_design.tsv`, `ruv_adjusted_protein_matrix.tsv`, `ruv_covariates_optimization.tsv`, `ruv_covariates_optimization_results.RDS`, `ruv_covariates_cancor.png`. |
| Transcriptomics | `trxRuvs` (`func_trx_ruv.R:140-157`) in `trx_ruvs` (`utils_workflow_trx_declarations.R:268-282`) | `W_1..W_k` in the design (:275); checkpoint `trx-ruvs`; product `trx_ruvs_science` with the display matrix `corrected`, W and k (:278); state `trx_ruvs`. File `ruv_normalised_log_values.tsv` (:557-558). |
| Transcriptomics | `trxRuviii` (`func_trx_ruv.R:165-178`) in `trx_ruviii_export` (:394-405) | `newY`, accepted under `snapshot$exports$trx_ruviii_export` (:418); no state. File `transcriptomics_normalized_data.tsv` (:643-649). |
| Lipid, metab | `ruvIII_C_Varying` (`func_lipid_s4_ruv_methods.R:397-935`, `func_metab_s4_ruv_methods.R:397-935`) in `omicsRuvOperation` (`utils_workflow_omics_declarations.R:305-329`) | Assay data replaced (:913); checkpoint `<prefix>-ruv-corrected`, optimisation product; state `lipid_ruv_corrected` / `metab_ruv_corrected`. No RUV-specific file. |
| Workbooks only | Protein-level calls in `DIA_workflow_starter.rmd:1345`, `DIA_workflow_experienced.rmd:867`, `DIA_workflow_limpa.rmd:1642`; metabolomics in `metabolomics_workflow_starter.rmd:580`, `metabolomics_workflow_experienced.rmd:677`; the peptide-level calls of section 3 | The objects of those workbooks. |

### 4.2 Readers

| Route | Reader | File and line | Values it shows today |
|---|---|---|---|
| Classic | QC figures `ruv_corrected_*.png` | `mod_prot_norm_qc_generation_helpers.R:228-281`, called at `utils_workflow_dia_products.R:646` | Corrected |
| Classic | Composite figure, column "RUV-Corrected" and the control plot | `utils_workflow_dia_norm.R:467-569`; called at `utils_workflow_dia_products.R:662` | Corrected; the control plot is of the data before RUV |
| Classic | Correlation filter | `diaCorrelationSource`, `utils_workflow_dia_products.R:1034-1039` | Corrected |
| Classic | `ruv_normalised_results_cln_with_replicates.tsv` | `utils_workflow_dia_products.R:1101-1108` | Corrected, after the correlation filter |
| Classic | Test input (limma) | `func_prot_da_model_methods.R:408-417` | Corrected, no W |
| Classic | Wide and long result tables, `log2norm.*` and `raw.*` columns | `func_prot_da_model_methods.R:535-574` | Corrected in both sets of columns (:565-566) |
| Classic | Glimma volcano's abundance panel; heatmap | `utils_workflow_dia_da.R:183-192, 201-251`; `func_prot_da_volcano_glimma.R:238-256`; `func_prot_da_heatmap.R:57` | Corrected |
| Classic | Result exports (`long_annot.tsv/.xlsx`, CSV download) | `func_prot_da_results_methods.R:58-140`; `utils_workflow_dia_da_views.R:335-356` | Statistics, with corrected per-sample columns |
| Classic | Results summary copy: `RUV_normalised_results.tsv`, figures | `func_general_filemgmt_results_summary_helpers.R:326-408, 643-680` | Corrected |
| Classic | Report | `inst/reports/proteomics/DIANN_report.rmd:222-289, 525-540, 646-674, 960-974` (LFQ and TMT reports alike) | Composite figure (corrected); RUV parameters as text |
| Classic | Enrichment | `utils_workflow_enrichment.R:845-866`; `func_general_s4_objects.R:165-193` | Result statistics only; no abundance matrix |
| Classic | Integration workbook | `inst/workbooks/integration/integration_workflow_starter.rmd:337` reads `ruv_normalized_results_cln_with_replicates.tsv` | A file name the app does not write (the app writes `ruv_normalised...`, with an s). Corrected only if the user supplies the file under that name. |
| Classic | Session export | `mod_prot_norm_session_helpers.R:19-101`; `utils_workflow_dia_norm.R:787-793`; `utils_workflow_dia_sessions.R:1607-1640` | RUV mode, k, the optimisation result, the `ruv_corrected_*` figures, the corrected table |
| Classic | Study parameters | `func_general_filemgmt_study_parameter_writer.R:99-145, 641`; `func_general_filemgmt_study_parameter_support.R:231-330` | Best percentage, k, scores, control count, metric, penalty, sample size, grouping, replicate structure. The RUV mode is not written in that section. |
| DIA_LIMPA | QC figures `ruv_corrected_*.png` | `utils_workflow_dia_products.R:988` | The adjusted matrix in both modes. In covariates mode this matrix is not what the test uses. |
| DIA_LIMPA | QC figures redrawn after a change of plot settings | `mod_prot_norm_qc_support_helpers.R:168-171` | Covariates mode: the uncorrected matrix, written under the `ruv_corrected_*` names. Correction mode: corrected. |
| DIA_LIMPA | Composite figure | None: `generateProtNormCompositeQcFigure` is called only on the classic route (`utils_workflow_dia_products.R:662`) | Not produced |
| DIA_LIMPA | `ruv_normalised_results_cln_with_replicates.tsv` | `utils_workflow_dia_products.R:1035-1036, 1105-1108` | Written by the correlation filter, which runs before RUV on this route: values from before RUV under the "ruv" name; copied to the results summary as `RUV_normalised_results.tsv` (`func_general_filemgmt_results_summary_helpers.R:381-386`) |
| DIA_LIMPA | Test input (`dpcDE`) | `func_prot_da_model_methods.R:366-406` | Covariates: uncorrected EList, W in the model. Correction: corrected EList, no W. |
| DIA_LIMPA | Result tables' per-sample columns, Glimma panel, heatmap, exports | As the classic rows | Covariates: uncorrected. Correction: corrected. |
| DIA_LIMPA | Report | `inst/reports/proteomics/DIANN_limpa_report.rmd:391-398, 528-554, 613-620, 894-918` | Text per mode (k, percentage, controls, replicate structure); the adjusted matrix labelled as QC only |
| DIA_LIMPA | Final object for the summary, study parameters and integration | `utils_workflow_dia_summary.R:84-92` | The first state present of `ruv_covariates`, `correlation_filtered`, `ruv_corrected`. In correction mode `correlation_filtered` comes first, so the object from before RUV is taken (order read here; the run was not made). |
| DIA_LIMPA | Study parameters | Writer :108-110, :641 | The optimisation result; a `[ruvCovariates]` section when the final state is `ruv_covariates` |
| DIA_LIMPA | Session archive | `utils_workflow_dia_sessions.R:1607-1640` | The state store; the `ruv_covariates_*` files and the adjusted matrix file are not in the archive |
| Trx | Stage figures, comparison grid, final PCA | `utils_workflow_trx_declarations.R:511-522, 565-595, 666-670` | The display matrix; the before-and-after figure shows the RUVIII export |
| Trx | `ruv_normalised_log_values.tsv` | :544-561 | The display matrix |
| Trx | Test input (edgeR) | :290-331, :806-810 | Original counts, W in the model |
| Trx | Result tables, volcanoes, summary text | :601-662; `func_trx_da.R:94-113` | Statistics only |
| Trx | Results summary and report | `func_general_filemgmt_results_summary_helpers.R:927-943`; `inst/reports/transcriptomics/transcriptomics_report.rmd:44-126` | The display matrix file and figures; W described as covariates |
| Trx | Integration export | :643-649; written to the integration directory (`utils_workflow_summary.R:119-127`) | The RUVIII matrix (levels shifted, section 1.4) |
| Trx | Integration workbook | `integration_workflow_starter.rmd:338` | That file, if the user puts it in the workbook's data directory |
| Trx | Study parameters and record | `trxRuvOutcomes`, :255-263; :440-453 | Mode, percentage range, k, percentage, control count, M source |
| Trx | Enrichment | None (`utils_omic_gui_config.R:2049`) | — |
| Lipid, metab | QC figures per assay, composite | `utils_workflow_omics_actions.R:340-369, 443-460`; `utils_workflow_omics_figures.R:188-230` | Corrected |
| Lipid, metab | Control plot | `mod_omic_norm_panels.R:430-432` | On the tab only; no file |
| Lipid, metab | Correlation filter, final PCA | `utils_workflow_omics_declarations.R:1049-1064` | Corrected |
| Lipid, metab | Test input (limma) | `func_lipid_da_model.R:251-283`; `func_metab_da_orchestration.R:7-46` | Corrected |
| Lipid, metab | `intensity.*` columns, Glimma, heatmaps, exports | `func_lipid_da_model.R:288-325`; `func_lipid_da_plotting.R:237-253, 441-475`; `func_metab_da_volcano_glimma.R:237-253`; `func_metab_da_heatmap.R:101-135`; `func_lipid_da_export.R:32-264` | Corrected |
| Lipid, metab | Results summary copy and report read `ruv_normalised_results.{tsv,RDS}` | `func_general_filemgmt_results_summary_helpers.R:495-523, 608-635`; `lipidomics_report.rmd:1091-1106`; `metabolomics_report.rmd:1050-1063` | No current code writes that file; the composite figure is corrected |
| Lipid, metab | Integration object `<omic>_<label>_final_s4.RDS` | `utils_omic_gui_config.R:1618-1621, 1835-1838`; `utils_workflow_omics_declarations.R:1390-1400` | Corrected |
| Metab | Integration workbook | `integration_workflow_starter.rmd:452` reads `workshop_data_final_metabolomics_s4_object.RDS` | Not the app's file name; lipidomics is not loaded by the workbook |
| Lipid, metab | Session export, study parameters | `utils_workflow_saved_exports.R:24-211`; `func_general_filemgmt_s4_helpers.R:146-189, 375-418` | The head state (corrected); per-assay optimisation; grouping, k and control count per assay |
| Lipid, metab | Enrichment | None (`utils_omic_gui_config.R:1589, 1817`) | — |

Recorded, not acted on (outside this ticket's scope): result tables carry `analysis_type = "RUV applied"` whatever the
mode (`func_prot_da_model_methods.R:486, 518, 546, 564`); `updateRuvParameters` stores the length of the control
vector, which is the number of proteins, as the number of controls (`func_prot_norm_optimization_helpers.R:298`), and
has no caller under `R/`.

## 5. Not read at the source, and errors in the ticket's table

**Not read at the source**

- Bradley 1978: the original. The criterion is given from three secondary papers that agree on the numbers.
- De Livera et al. 2012, Redestig et al. 2009 (CCMN), Han & Li 2022: no open copy.
- Supplements: Poulos 2020's Supplementary Notes (which tests reproduced the plasma study), Molania 2019's
  Supplementary Method (the RUV-III algebra, so whether that paper's formula centres the controls is unknown),
  Molania 2023's supplement, Gagnon-Bartsch & Speed 2012's supplement.
- Li et al. 2023: read in pieces; the middle step of ComBat+Cor only as a paraphrase.
- `RUV4` and `RUVinv` help: argument and details sections only. `RUVinv`'s arguments not inspected.
- Morris 2019 Table 6: the standard-error formula for a rejection rate was reconstructed from the worked examples.
- `ruv` was read on the CRAN mirror at 0.9.7.2; the pinned library holds 0.9.7.1. `RUVIII`, `RUV1` and
  `ruv_cancorplot` were printed from the installed copy and match what is quoted.
- Not run: no package route was executed. The level statements for lipid, metab and the transcriptomics export rest on
  the code and on the toy check of the two reference functions.
- `AGENTS.md` requires four files under `conventions/` to be read; the directory is absent on `kit-123`.

**Errors and gaps in the ticket's table ("The five implementations today")**

1. *Reference, controls: "RUV-III-C paper: spike-ins".* Wrong. The paper used peptides from yeast proteins rescaled to
   equal means, and in its second dataset 8,055 human peptides detected in every sample; it says endogenous features
   are generally the better choice than spike-ins.
2. *Reference, level: "Kept".* Wrong for the two packages W1 compares against. `ruv::RUVIII` and
   `RUVIIIC::RUVIII_C_Varying` do not centre and shift each feature's level; only the PRPS authors' code keeps it.
   Consequence for W1: on kit-123 the classic method no longer equals `RUVIII_C_Varying` and `computeRuvCovariates` no
   longer equals `ruv::RUVIII`'s `newY`; each differs by a constant per feature. "Results equal" has to be stated up
   to that constant, or against the PRPS formula.
3. *Reference, test: "Patterns in the model; corrected values for display (RUVSeq, limma)".* One side only. The RUV-III
   paper (Molania 2019) tests corrected values with limma and a design of the biological factors, which is what the
   classic, lipid and metab routes do.
4. *Reference, k: "No rule; small, by diagnostics".* "No rule" is right. "Small" holds for RUVSeq's examples (k = 1)
   and not for the RUV-III family (20 and 11 in the RUV-III-C paper; 6, 10 and the maximum in Molania 2019; 20 and 12
   in the PRPS code; 5 in hRUV).
5. *Reference, replicates: "Technical replicates".* True of the RUV-III-C and RUV-III papers. The RUVs vignette uses
   the treatment groups, and PRPS uses pseudo-replicates within a biological group.
6. *"Every omic also chooses its controls and its k by the same objective: the largest gap".* Not for k. Classic,
   lipid, metab and the peptide method take the elbow (the smallest K within 5% of the largest gap); DIA_LIMPA and
   transcriptomics take the K at the largest gap. The percentage is chosen by the gap times a penalty on k. The
   q-value method also differs: `qvalue` in the classic automatic sweep, BH elsewhere.
7. *k, DIA_LIMPA "Bounded, 2 to 5"; transcriptomics "Bounded, 3 or 4".* These are the caps. k runs from 1 to the cap
   (and below the number of controls for DIA_LIMPA).
8. *Level, transcriptomics: "Kept".* True of the display matrix. The RUVIII export, which is the matrix written for
   integration, is `ruv::RUVIII`'s output and is shifted. The code row also omits `ruv::RUVIII` for this omic.
9. *Level, lipid and metab: "Same call as classic".* No longer the same behaviour on kit-123: R2 restores levels in the
   proteomics method only. Lipid, metab and the peptide method return what `RUVIII_C_Varying` returns.
10. *k, classic: "(KIT-123 R1 bounds it)".* Right for the automatic sweep. Manual mode has no upper limit, and the
    lipid, metab and peptide sweeps are not bounded by R1.
11. *Scope.* The table has no column for the peptide-level method; it is reachable only from three workbooks.
12. *Standards.* "Measures are fixed in advance" is attributed to Weber et al. 2019 and is not in that paper. Bradley's
    band is two-sided and W3 rule 1 applies only its upper limit.
