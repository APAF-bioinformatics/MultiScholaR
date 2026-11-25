# KIT-267 wording (drafts for the owner's wording pass, VAL-41 and VAL-311)

Texts that changed or are new. The wording is the owner's; the tests require only the mechanism.

| Where | Text (draft) |
|---|---|
| Wording key `mark_organism` (after the labels of STRING species, KEGG organism, Reactome organism) | human; change for other species |
| `settled_by` of the eight (script users of `integrationSettings()$enrichment`) | "KIT-141 (chosen 2026-11-21 under INT-18; decided 2026-11-26, VAL-310; methods ledger, the ranked list)"; "... under INT-19; decided 2026-11-26, VAL-310; methods ledger, KEGG and Reactome rules)"; STRING species: "KIT-141 (human, 9606, as the starting value; decided 2026-11-26, VAL-310)" |
| NEWS | The Enrichment tab's defaults are confirmed decisions; the three organism settings say they default to human and can be changed for other species (KIT-240, KIT-267). |
| Testers' guide, section 6 | the sentence about eight Enrichment settings marked "default not yet confirmed" is removed |

## Items for the owner's wording pass (none blocks)

1. The double parenthesis after the three organism labels, e.g. "KEGG organism code (hsa is human, mmu mouse) (human; change for other species)": shorten the three labels so the mark carries the hint alone, or use a help line under each input.
2. Whether the other Enrichment settings should show "(decided)" as the Model and Results tabs do (recommendation: no mark).
3. The integration report's status legend still prints the row "default not yet confirmed" (`inst/reports/integration/integration_report.rmd:51-52`; `R/utils_integration_record.R:38`, `:81`); after this ticket no Model, Results or Enrichment setting can carry that status.
