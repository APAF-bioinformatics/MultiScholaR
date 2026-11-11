# KIT-116 W2, A101-10: the mutation check on a proteomics lane (2026-11-11)

KIT-101's check (b) was run on the trx lane only. Here it is run on `prot_dia_limpa`, KIT-101's procedure: a scratch copy of the tree at `3c2f9866` (`git archive`), one default changed, the lane alone against its committed oracle (`tests/testdata/gui-entry/prot_dia_limpa.json`, frozen at `08e9b25e`), `load_all`.

| Attempt | Change in the copy | Result |
|---|---|---|
| 1 | `R/mod_omic_da_panels.R:190, :522`, the DA heatmap's `top_n` default 50 to 40 | **Not detected:** 38/0. The lane does not capture the heatmap view: the proteomics DA views are covered by `testServer` blocks only (audit A106-03, F106-01; KIT-116's e2e item). A valid result, kept. |
| 2 (the item's one correction pass) | `R/mod_prot_qc_protein_normalisation.R:47`, the DIA_LIMPA normalisation default `none` to `cyclicloess` | **Detected:** 15 pass, 23 fail. |

**Attempt 2, by step:**
- From `import` to `qc_protein_dpc_quant`: only `dynamic_qc_tabs` differs (the rendered tab shows the default).
- At `qc_protein_normalisation`: the first data difference, the `protein_normalised` state.
- After it: every later state (`peptide_count_filtered`, `correlation_filtered`, the RUV states), the RUV artifacts and table, the DA tables and figures, and the summary's files.

**Read:** the lane detects a changed scientific default, at the step that uses it and in everything downstream. It does not detect a change confined to the DA heatmap view. That gap is A106-03's, and KIT-116's DA-views e2e lane is meant to close it.
