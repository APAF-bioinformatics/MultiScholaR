# KIT-173 results, in plain words

Gene Ontology parent terms on the proteomics enrichment route for organisms outside the 13 that g:Profiler serves.
Figures and runs are in `outcome.md`, `runs.tsv` and the tables named below.

## 1. Summary

- **Agreement with g:Profiler (T5).** On the reviewed yeast proteome (UniProt 2026_04), the route's propagated sets agree
  with g:Profiler's own propagated sets over 2,703 terms: median containment 0.933 (interquartile range 0.833 to 1.000)
  and median precision 1.000. Without propagation the median containment is 0.073. All three conditions of the rule
  hold (`t5-summary.tsv`). With is a and part of only, containment is 0.893 over all terms, 0.935 over the 1,592 terms
  that no regulates edge reaches and 0.810 over the other 1,111 (`t5-two-relation.tsv`).
- **Effect (9.1).** On one proteome (reviewed E. coli K-12, 4,024 annotated proteins) with synthetic lists, propagation
  makes about four times as many terms testable (314 against 1,249). Planted process terms are found in every list
  with parent terms. Without them, two of the three planted terms cannot be tested at all, because they have too few
  directly annotated proteins. On random lists, without parent terms, one list had 120 terms significant by q-value
  and none by Benjamini-Hochberg; with parent terms the two flags never differed by more than 3 terms. The figures
  describe this proteome and these lists, not power or a false-discovery rate (`outcome.md`, Effect).
- **Replay (A4, as restated under VAL-61).** Inside a `geneID` cell, the order of the proteins was set by the R process
  (enrichit's hash set over `Rcpp::unique()`), on the base as well. The route now lists them alphabetically, and records
  made after the fix replay byte-identically across processes. Records made before the fix replay with every value
  equal and their `geneID` cells equal as sets.

## 2. T5's thresholds

The rule (median containment and precision at least 0.90, and containment at least 0.20 above that without
propagation) was set on this test's own fixture, before the pilot. It is a gross-error bound: it would catch a build
that does not propagate, or one that propagates across aspects or to every term. It is not an estimate of accuracy.
g:Profiler builds its sets independently, from its own annotation import (Ensembl BioMart, data e114_eg62_p19) and its
own handling of relations. UniProt, GO.db and g:Profiler releases also differ by a few months. Release drift was
small: no g:Profiler term is missing from GO.db 3.23.1, and no g:Profiler gene was lost in mapping.

## 3. The lowest-agreement terms (`t5-lowest.tsv`)

The table lists the 20 terms with the lowest containment and the 20 with the lowest precision, under the shipped
value. Causes as far as the records show. These are inferred from the term names and sizes; no annotation-level
comparison was made.

- **Molecular-function kinase terms with about 92 g:Profiler genes against 1 to 3 here** (GO:0004679 AMP-activated
  protein kinase activity, 1 against 99; GO:0004694 eukaryotic translation initiation factor 2alpha kinase activity,
  1 against 92; GO:0044025, GO:0140998, GO:0072354, GO:0140996, GO:0035173, histone kinase activities, 92 each). The
  shared size of 92 suggests that g:Profiler's annotation source assigns a whole family (yeast protein kinases) to
  each specific kinase term. That would be an annotation-source difference between Ensembl's import and UniProt's
  `go_id`, not a relation difference: precision is 1.000, so every protein here is also in g:Profiler's set. Not
  verified.
- **Regulation-of-activity process terms** (GO:0071900 regulation of protein serine/threonine kinase activity,
  GO:0000079, GO:0051336, GO:0043087, GO:0051345, GO:0043547; containment 0.03 to 0.06, precision 0.5 to 0.67).
  g:Profiler's larger sets are consistent with propagation from molecular-function regulator terms into these
  process terms. Those links cross aspects and are absent from go-basic.obo, which GO.db is built from, so they are
  not followed here. This is a relation-source difference. Not verified.
- **Small salvage and metabolic terms with 0 overlap** (GO:0043102 amino acid salvage and GO:0071267 L-methionine
  salvage, 1 protein here against 10, none shared). Different proteins are annotated in the two sources: an
  annotation or release difference.
- **Low-precision process terms** (GO:0071265, GO:0000727, GO:0009083, GO:0000097 and others; precision 0.27 to 0.63).
  Here the route's set is larger than g:Profiler's. That is consistent with UniProt's `go_id` carrying annotations
  that g:Profiler's import does not. GO:0000727 has containment 1.0 (63 here against 26), so the extra proteins come
  from direct annotations or their descendants in UniProt. Not verified.
- **Other terms in the list**: no cause identifiable from the records.

## 4. Methods sentence (proposed, VAL-41; the owner edits)

"For organisms outside g:Profiler's set, Gene Ontology protein sets were built from UniProt annotations and
propagated to parent terms over the is a, part of, regulates, positively regulates and negatively regulates relations
of GO.db 3.23.1 (go-basic.obo of 2026-01-23), so that a biological-process term also holds proteins annotated as
regulating that process (is a and part of only is available as an option); the three top-level terms were not
tested. Over-representation was tested with clusterProfiler against the measured proteins with at least one
annotation below the top-level terms, over terms with 10 to 500 such proteins and at least one protein in the list,
the three aspects together as one family per comparison and direction; a term was called significant when its Storey
q-value was below 0.05, with Benjamini-Hochberg adjusted p-values reported beside it. GO terms are nested, so the
tests are positively dependent and a significant term's parents tend to be significant too: significant terms are
reported as related groups, not as independent findings. Storey's q-value assumes nearly independent tests; on 20
random lists from the E. coli proteome 1 list in each setting had any significant term (95 % interval 0.1 % to 25 %),
and without parent terms one random list had 120 terms significant by q-value and none by Benjamini-Hochberg."
