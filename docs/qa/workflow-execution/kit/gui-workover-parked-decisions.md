# GUI workover: decisions parked for the user

The user, 2026-11-07: "continue with this series agentically until the GUI is done. any fixes along the way - do it as long as the scientific methodology doesnt change. anything that needs my input park for end of the chain."

Each entry gives: the ticket, the question, the evidence (links), the options with a recommendation, and what stays blocked until it is answered. Nothing parked here has been acted on.

| # | Ticket | Question | Evidence | Options (recommended first) | Blocks |
|---|---|---|---|---|---|
| P-01 | KIT-114 | The proteomics DA q-value method (`qvalue()` against BH at small n) and RUV-III correction mode at small n (88 false positives at n = 8). Both are methodology. | `kit-101/diagnosis/f19.txt`; KIT-114 (a), (b) | Measure on synthetic truth across sizes, then decide | KIT-114 only |
| P-02 | KIT-114 | DIA_LIMPA's normalisation point: E-level (F-23, approximate for quantile and scale) or precursor level before DPC (limpa's documented point) | `kit-101/findings.md` F-23; KIT-114 (c) | Keep E-level, default `none` (today); measure the precursor level | KIT-114 only |
| P-03 | KIT-114 | RUV negative-control selection with missing values: restrict to complete proteins, or impute for the k count | `kit-101/findings.md` F-27; KIT-114 (d) | Message now (being done); measure the method | KIT-114 only |
| P-04 | real-data gate | KIT-054 `kc_lfq` and `kc_lfq_dotted` oracles were recorded before F-18 (a 0 kept as log2 0). Re-run at the final real-data gate. | `kit-101/findings.md` F-18 | Re-record at the gate (as decided) | nothing in the GUI series |
| P-05 | KIT-102 | Non-numeric text in a proteomics quantity column (Spectronaut "Filtered", decimal-comma exports, Excel error cells). DIA-NN, MaxQuant and FragPipe silently read it as NA (`as.numeric`, `func_prot_import_readers.R:43,203,330`). After KIT-102, Spectronaut and PD-TMT convert only numeric-like text, and anything else fails at the design step as today. Lipid, metab and trx refuse the import with a message. | KIT-102 parity section; `kit-102/parity.md` | Keep each reader as it is (no scientific change) / refuse at import with a message naming the column, for every reader / map known vendor tokens to NA with a notification | nothing in the GUI series |
| P-06 | KIT-102 | Trx's `column_mapping` key is `lipid_id_col`, a misnomer kept because `trxDesignObject` and the trx98a oracles read it. | KIT-102 | Keep it (a journey and state contract) / rename it to `gene_id_col` and re-record trx98a | nothing |
