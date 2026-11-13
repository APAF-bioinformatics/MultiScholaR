# KIT-114: measurements (a), (e) and (f) (2026-11-13)

**Outcome:** (a), (e) and (f) measured; no default changed; each decision is parked for the user (`parked-20261113.md` rows 6 to 8). (b), (c) and (d) are deferred, each needing its own simulator design (`tickets/kit/KIT-114.md`, revision 1.2). Parameters and rules were fixed in the ticket before any run. Scripts: `sim/`; per-seed tables `*.tsv`, summaries `*-summary.tsv`.

## (a) `qvalue()` against Benjamini-Hochberg

Through the proteomics DA kernel (`runTestsContrasts`), whose significance column is `fdr_qvalue < 0.05`. 2,000 centred proteins, n = 3, 4, 6, 10 per group, 5% or 20% truly different, |log2FC| 0.5 or 1.0; 16 cells × 20 seeds (`a_qvalue_bh-summary.tsv`).

- **Both control the FDR in every cell** (achieved ≤ 0.05 + 2 SE). Where calls are numerous, qvalue sits at about 0.05 and BH at about 0.04 (BH is conservative by the null share).
- **qvalue's power gain** is small to moderate: +0 to 3 points at 5% different; +2 to 5 points at 20% (largest at n = 4, |log2FC| 0.5: 8.2% against 3.4%).
- In cells with almost no power (n = 3, |log2FC| 0.5) a few calls make the false-discovery proportion jumpy (0.10 to 0.15 ± 0.07 to 0.08); both columns behave the same there.
- **Rule outcome:** both control everywhere, so the recommendation is to **keep `qvalue()` as the significance default and report both** (both columns are already in the results table).
- **Caveats:** a centred simulator with limma's own variance model; it does not test RUV, imputation or the q-value's π0 estimate under correlated proteins. The audit's earlier "88 false positives" came from an uncentred, single-seed fixture and is superseded.

## (e) Lipid and metab zeros: kept (log2(x + 1), today) against missing

Through the lipid DA kernel (`runTestsContrastsLipidDA`; the metab kernel is the same code). 500 features, n = 3, 6, 10, 20% different with |log2FC| 1.0, left-censoring written as 0 for 10% or 30% of values; both treatments on the same features (at least 2 non-zero values per group, as the QC filter leaves them; a correction made after the first run failed on that point, before any result); 6 cells × 20 seeds (`e_zeros-summary.tsv`).

| n | zeros | FDR kept / missing | power kept / missing | mean \|log2FC error\| as share of the true effect, kept / missing |
|---|---|---|---|---|
| 3 | 10% | 0.032 / 0.053 | 0.64 / 0.78 | 81% / 20% |
| 3 | 30% | 0.021 / 0.050 | 0.55 / 0.75 | 142% / 20% |
| 6 | 10% | 0.022 / 0.044 | 0.72 / 0.93 | 105% / 15% |
| 6 | 30% | 0.015 / 0.049 | 0.48 / 0.91 | 223% / 16% |
| 10 | 10% | 0.026 / 0.051 | 0.69 / 0.97 | 116% / 12% |
| 10 | 30% | 0.023 / 0.054 | 0.46 / 0.95 | 244% / 13% |

- Both control the FDR (keeping zeros is conservative).
- **Keeping zeros loses a large part of the power** (14 to 49 points) and **distorts individual fold changes**: the mean absolute error is 81% to 244% of the true effect against 12% to 20% with zeros as missing. The *signed* mean bias is near 0 for both (the errors cancel across features).
- **Rule outcome:** the rule said "biases log2FC by more than 25% of the true effect". Read as per-feature error it is met by a wide margin; read as the signed mean bias it is not. Reported both ways. **Recommendation: read 0 as missing for lipid and metab**, for power and estimate accuracy, matching proteomics (F-18) and the DIA-NN zero rule (KIT-118). The same decision should settle negatives (lipid clamps them to 0 before `log2(x + offset)`).
- **Caveats:** no normalisation or ITSD step in the simulation (centred data); real normalisation would also be distorted by zeros, which this does not measure.

## (f) Enrichment: id stripping and organism matching

- **Id stripping** (`func_general_enrichment_orchestration.R:121` strips `:.*` from the query ids; the background at `:217` does not). Enrichment uses `uniprot_acc`. Standard forms (`P12345`, `P12345-2`, `P12345;Q67890`, `sp|P12345|ALBU_HUMAN`) contain no colon and pass unchanged, and nothing in the package builds an id with a colon. **No result changes for standard accessions.** A colon-bearing id (a custom FASTA header) would leave the query and fall out of the background match. **Recommendation:** strip both the same way, in KIT-119. It changes results only for such ids.
- **Organism matching.** Proteomics' lookup is an exact taxon-id match (`utils_workflow_enrichment.R:220-240`, `utils_workflow_dia_enrichment.R:384-404`): no issue. The substring match is in `runReactomeEnrichment` (`func_multiomics_enrich_reactome_helpers.R:72`, exported; used by the integration notebook, not the app). An input that matches several Reactome species takes the first (for example "saccharomyces" matches both *S. cerevisiae* and *Schizosaccharomyces pombe*), and no match falls back silently to *Homo sapiens*. **Recommendation:** exact match, with an error on an ambiguous or unknown organism. This changes results only where the input was ambiguous or unknown.
