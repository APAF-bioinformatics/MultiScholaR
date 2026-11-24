# KIT-205: user-facing texts changed (drafts for the owner's wording pass, VAL-41; plan section 5)

Every text below is proposed wording. Old and new are given in full. Nothing else on the tab changed.

## The tab (KIT-203's texts, SD-1 (b): "measured" becomes "tested")

| Where | Old | New |
|---|---|---|
| Significance rule selector, default choice (`ENRICHMENT_SIGNIFICANCE_RULES$label`; also printed in the completed status's "Significance:" line and in the download summary) | Benjamini-Hochberg over every term of 10 to 500 measured proteins (default) | Benjamini-Hochberg over every term of 10 to 500 tested proteins (default) |
| Default rule's note beside the selector (`$note`), first sentence | Each term's p-value is adjusted together with every Gene Ontology term (biological process, molecular function and cellular component together) that has 10 to 500 of your measured proteins, including terms with none of the listed proteins (they count as p = 1). | Each term's p-value is adjusted together with every Gene Ontology term (biological process, molecular function and cellular component together) that has 10 to 500 of your tested proteins, including terms with none of the listed proteins (they count as p = 1). |
| Run notice, clusterProfiler route, default rule (`$notice`) | For this organism, terms are called significant by the Benjamini-Hochberg adjusted p-value over every term of 10 to 500 measured proteins. | For this organism, terms are called significant by the Benjamini-Hochberg adjusted p-value over every term of 10 to 500 tested proteins. |
| clusterProfiler results table legend (`R/mod_omic_enrich.R`), second sentence | p_adjust_family: Benjamini-Hochberg over every term of 10 to 500 measured proteins. | p_adjust_family: Benjamini-Hochberg over every term of 10 to 500 tested proteins. |

## The four proteomics report templates (DIANN, DIANN_limpa, LFQ, TMT)

"With the rule" means a record whose enrichment step carries `enrichment_background = "tested"` (every run from this
change on), or a report rendered without a record. "Earlier record" means a record without it (made before this change,
with every row of the table as the background).

| Sentence | Old (all cases) | New, with the rule | New, earlier record |
|---|---|---|---|
| clusterProfiler family, default rule | ...: terms with 10 to 500 annotated measured proteins, whether or not a listed protein falls in them (terms without one enter with p = 1). | ...: terms with 10 to 500 annotated tested proteins, whether or not a listed protein falls in them (terms without one enter with p = 1). | ...: terms with 10 to 500 annotated proteins of that comparison's differential-abundance table, whether or not a listed protein falls in them (terms without one enter with p = 1). |
| clusterProfiler family, q-value rule | ...: terms with 10 to 500 annotated measured proteins and at least one protein in the list. | ...: terms with 10 to 500 annotated tested proteins and at least one protein in the list. | ...: terms with 10 to 500 annotated proteins of that comparison's differential-abundance table and at least one protein in the list. |
| clusterProfiler background | The background was the measured proteins with at least one Gene Ontology annotation{background_clause}. | The background was the proteins tested in that comparison (those with a p-value) with at least one Gene Ontology annotation{background_clause}. | The background was every protein listed in that comparison's differential-abundance table, including any without a p-value, with at least one Gene Ontology annotation{background_clause}. |
| g:Profiler background bullet (with KIT-208's clause, E-7) | • The background was the identifiers of each comparison's differential-abundance table (gene names where the table has them, otherwise protein accessions) that g:Profiler recognises (a custom background). | • The background was the identifiers of the proteins tested in each comparison (those with a p-value in its differential-abundance table; the first gene name of each protein in the protein annotation; proteins without a gene name were left out of the lists and the background) that g:Profiler recognises (a custom background). | • The background was the identifiers of each comparison's differential-abundance table (the first gene name of each protein in the protein annotation; proteins without a gene name were left out of the lists and the background) that g:Profiler recognises (a custom background). |

`{background_clause}` is KIT-173's " (on a run with parent terms, at least one below the three top-level terms)" on a
propagated run, else empty (unchanged).

## Help, NEWS and the log

- `processEnrichments` help, `significance_rule`: "over every term of 10 to 500 measured proteins" becomes "over every
  term of 10 to 500 tested proteins". New Details paragraph: "The background (universe) of each comparison is the
  proteins tested in that comparison (KIT-205, VAL-129 (2)): before the up and down lists and the background are formed,
  the comparison's table is reduced to its rows with a finite `fdr_qvalue` and, where the table has a `raw_pvalue`
  column, a finite `raw_pvalue` (a table without `raw_pvalue` is reduced by `fdr_qvalue` alone). A protein without a
  p-value could never be in a list, so it is not counted in the background either. One background per comparison
  serves both lists, on both routes; the number of rows left out is written to the log. Analyses recorded before this
  rule replay with every row."
- NEWS, new entry under Scientific changes: the ticket's text with K205-Q03's worked example (see `NEWS.md`).
- NEWS, KIT-203's entry: "over every term of 10 to 500 measured proteins" becomes "tested proteins".
- Log line (both routes; no notice, ED-5): old "Using %d unique proteins as background for enrichment analysis"
  (g:Profiler only); new "Using %d unique proteins as background for enrichment analysis (%d rows without a p-value
  left out)".
- The declaration's `global_statistics` text (validated, not shown or recorded): old "the background universe and the
  significant set are the complete DA table of each contrast"; new "the background universe is the rows of each
  contrast's DA table with a p-value; the lists are drawn from it".
