# KIT-201 wording drafts for the owner (VAL-41)

Every public text below is a draft for the owner's wording pass; none is final.

1. **NEWS** (`NEWS.md`, the KIT-201 entry under "Scientific changes"): as committed, with AUDIT-KIT-201-impl I02
   (lipid/metab figures are features whose true change equals the threshold), I03 (the two measured facts: proteomics
   at threshold 0 with abundance-dependent missing values stays above target; without an abundance trend the new test
   is slightly more liberal) and I04 (lipidomics and metabolomics are one experiment). The audit notes the owner may add
   "(shown when the RUV design panel is enabled)" to the Legacy-row sentence.
2. **Ledger** (`docs/methods-decision-ledger.md`, the two `VAL-<merge>` rows), including DAS-015's new public wording:
   "limma linear models with empirical Bayes moderation toward an abundance-dependent prior with robust hyperparameters
   (trend and robust), the same moderation in TREAT; TREAT at a log fold-change threshold of 0 by default, which at 0
   equals the moderated t-test; significance at q < 0.05".
3. **In-app Legacy row** (`R/utils_omic_gui_config.R`, classic proteomics): "For comparison with the approach of earlier
   versions; analyses saved with earlier versions replay unchanged. False discoveries ranged from 3% to 39% across the
   scenarios tested, against a target of 5%; 10.1% in the reference scenario."
4. **Methods** (ticket section 11; Manuscripts repository, not edited here): "Proteomics, metabolomics and lipidomics use
   limma linear models with empirical Bayes moderation of the variances toward an abundance-dependent prior
   (limma-trend; Law et al. 2014, Genome Biol 15:R29) with robust hyperparameter estimation (Phipson et al. 2016, Ann
   Appl Stat 10:946-963); fold-change thresholds use TREAT (McCarthy and Smyth 2009, Bioinformatics 25:765-771) with
   the same moderation." Limitation to add (I03): with 20 % abundance-dependent missing values the proteomics test at a
   threshold of 0 remained slightly liberal for low-abundance proteins in simulation.
5. **KIT-177's Methods limitation, corrected (AUDIT-KIT-201-impl I11; draft, Manuscripts repository, not edited here).**
   Replace "each protein's variance is moderated toward one prior variance shared by all proteins; no
   abundance-dependent trend is fitted, so the larger variance of log values at low abundance is not modelled" with:

   > each protein's variance is moderated toward a prior that follows its average log abundance, with robust estimation
   > of the prior (limma-trend; Law et al. 2014; Phipson et al. 2016). In simulation with 20 % abundance-dependent
   > missing values the test at a log fold-change threshold of 0 remained slightly liberal for low-abundance proteins
   > (null rate up to 5.8 %, false discovery proportion up to 5.7 % at a 5 % target; KIT-201 M1). R0's figures [9.1]
   > were measured with the earlier test, which used one prior variance for all proteins.
6. **Reports re-rendered from analyses saved before KIT-201** (audit section 8, choice (a) recommended, as implemented):
   the NEWS disclosure only; (b) would be a template change reading the record.
