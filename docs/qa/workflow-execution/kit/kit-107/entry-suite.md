# KIT-107 entry-point suite

## D6 test 2: the F-16 lane check fails first

The proteomics entry helper's `prot_dia` enrichment step gains one expectation: the analysis record built from the owner lists no incomplete `enrichment_completed`. On the start tree (`8e7294b6`, `git archive`, with only the new helper copied in), `prot_dia` run alone fails exactly that expectation ("KIT-107 F-16: the analysis record's enrichment step is complete"); every other expectation passes. Log: [lifecycle/start-prot_dia.txt](lifecycle/start-prot_dia.txt).

## Simulation Setup (run 2): prot_dia alone on the merged tree `2b1e02c8`

The lane was run through `kit-105/lifecycle/lifecycle.R` as committed, with strict skips and `R_LIBS` set to the KIT-106 corrected build (J57eHO). Log: [lifecycle/merged-prot_dia.txt](lifecycle/merged-prot_dia.txt).

**Result:** identical to its oracle (the ListReporter flag is TRUE) in 241 s, and F-16's lane check passes.

**Traced calls:**

| Function | Calls |
|---|---|
| `mod_omic_enrich_server` | 1 |
| `runWorkflowEnrichment` | 1 |
| `runDIAReferenceEnrichment` | 1 |
| `saveWorkflowEnrichmentState` | 1 |
| `resolveProtEnrichUniprotAnnotations` | 1 |
| `buildProtEnrichAllContrastResults` | 1 |
| `restoreWorkflowEnrichment` | 1 |

The ticket expected `restoreWorkflowEnrichment` to be called 0 times. It is called once because of F107-02's fix: the selected-contrast observer asks the core once before the Run, the core returns NULL (no state yet), and nothing is restored. The ticket's expectation predates the fix; the lane is identical either way. W3 measured the same.
