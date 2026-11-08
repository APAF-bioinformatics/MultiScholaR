# KIT-105 findings

Append only. Ids are F105-NN. A fix that does not change the methodology is made here, with a test that fails first. Anything else is parked.

| Id | Where | What | Status |
|---|---|---|---|
| F105-01 | `R/mod_lipid_norm_support_helpers.R:4-58, :90-140` (lipid and trx normalisation tabs) | `registerLipidNormStaticQcImageOutputs` and `registerLipidNormAssayLabelOutputs` pass lazy promises from a for loop into the renderer closures, so they are evaluated after the loop ends. Every lipid QC image slot shows the last spec (assay 2, ruv_corrected, correlation), and every label shows assay 2. Trx registers the same functions, so its normalisation tab shows no figures and its labels read "Assay 2: (detecting...)". Display only; no result changes. The frozen lipid_canonical and trx_screen oracles contain the defect. | **fixing (lead decision under the standing authority):** `force()` the renderer arguments, as its own commit. The gated W1 block goes green; lipid_canonical and trx_screen are re-frozen once, with attribution in [lane-changes.md](lane-changes.md). |
| F105-02 | skip-correlation's store commit (lipid `handleLipidNormSkipCorrelationFilter`, metab `dispatchMetabNormSkipCorrelation`) | An error in the commit is unhandled in both omics, so shiny ends the user's session. No result changes. W1's (c) capture freezes today's behaviour. | **assigned to W4 (the shared frame):** catch it, notify, keep the state. The one (c) capture is re-frozen with attribution. |
