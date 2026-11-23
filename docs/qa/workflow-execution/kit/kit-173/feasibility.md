# KIT-173 feasibility (F1, F2)

Ticket revision 1.5. F1 ran on `run-tests` (R 4.6.1, reference BLAS, ubuntu-22.04) as R0 (run 45119655641, commit
36a34cab) and its correction pass R0b (commit 32ffcf04); the run ids, minutes and core-hours are in `runs.tsv`. The
script is `f1-feasibility.R`; its outputs are under `f1/` (R0b; R0's where only R0 has them, named `r0-`). F2 is a
reading of the code on the branch and of g:Profiler's published documentation as fetched by F1 (no call from the
laptop).

## F1: observed values

__F1_TABLE__

## F1 stop rules

None triggered (`f1/stop-rules.tsv`, R0 and R0b):

- UniProt's `go_id` holds no ancestors: for P00330, P00549 and P02829 every UniProt GO id is among QuickGO's direct
  (non-NOT) annotations of the entry (13/13, 9/9, 25/25); none of the 72, 103 and 107 ancestors of those terms that
  are not themselves annotated appears in `go_id`.
- UniProt's `go_id` holds no NOT-qualified term: the three QuickGO `NOT|involved_in` annotations checked (P27697
  GO:0006468, P34227 GO:0050821, P38013 GO:0050821; none with a positive annotation to the same term) are absent
  from the entry's `go_id`.
- Every value of section 2.2 was reproduced (edge counts per label and aspect; pairs 464,934 / 352,532 BP, 64,418 MF,
  42,075 CC; `GO*ANCESTOR` identical to the five-relation closure in each aspect; no edge to another aspect except
  each root's `isa` edge to `all`; GO:0006110 18 against 75 ancestors; 3,297 secondary and 9,861 obsolete ids).
- No fetch failed (every fetch answered 200 at its first attempt in R0 and R0b).

## F1: what R0 found that the ticket did not expect, and what was done

1. **No root-only protein exists in either fixture.** Of the reviewed entries, none has GO ids that are all aspect
   roots (yeast 0 of 6,052 with GO, E. coli 0 of 4,024), and none has GO:0008150 alone. R0's script stopped at the
   T7 selection on this. T7 and T8 need exactly one such protein (ticket 10: "T7's fixture includes exactly one,
   chosen in F1"). R0b makes it synthetic and says so: an E. coli entry without any GO annotation, drawn with seed
   20261123 + 8, is given GO:0008150 alone in the T7 fixture (as a curator's "no biological data available"
   annotation would be). The README names the entry and the edit. Consequence for 5.1 and 9.1: on these two reviewed
   proteomes `go_root_only_proteins` is 0 for real data; the count matters for unreviewed entries, which carry such
   annotations, and is still recorded per run.
2. **The GMT's term names.** R0's upload scan withheld the GO GMT (a term name matched one of the scan's patterns).
   R0b commits the GMT with each line's description field replaced by the term id (ids and genes only; T5 reads only
   those). The SHA-256 of the file as g:Profiler serves it and of the committed file are both recorded.

## F2: answers

| Question | Answer | Where |
|---|---|---|
| Carrier of the gene-set attributes to the step outcomes without changing `EnrichmentResults` | None in the result object. The stop rule's fallback is used: `runWorkflowEnrichment` computes the outcomes (`workflowEnrichmentGoOutcomes`) from the annotation table it passes to the kernel and from `buildUniprotGoGeneSets()` called on it, the same call `processEnrichments` makes. Cost: the route's gene-set build runs twice on a propagated run (T6 measures one build) | `R/utils_workflow_enrichment.R` |
| Does `getUniprotAnnotations`' attribute reach `processEnrichments` on the application route? | Yes. The table the resolver returns (bus value, `<source_dir>/uniprot_dat_cln.RDS` via `saveRDS`/`readRDS`, or `getUniprotAnnotations`, whose cache is an RDS) is passed unchanged as `go_annotations` (`runWorkflowEnrichment` -> `buildProtEnrichProcessEnrichmentsArgs` -> the operation's kernel). A table cached before this change has no attribute and records `"not recorded"` | `R/utils_workflow_dia_enrichment.R` (`resolveProtEnrichUniprotAnnotations`), `R/utils_workflow_enrichment.R`, `R/utils_workflow_dia_da.R` |
| Through `getUniprotAnnotationsFull`? | Only the Ensembl branch calls `merge()`, which drops attributes; the attribute is re-attached after it. The other branch only adds attributes | `R/func_prot_annotation_uniprot_full_helpers.R` |
| Does `reviewed` reach the route? | Yes: `directUniprotDownload` asks for `reviewed`; `read.delim` names the column `Reviewed` (values `reviewed`, `unreviewed`); `uniprotGoIdToTerm` and `standardizeUniprotColumns` keep it. `go_unreviewed_share` is recorded | `R/func_prot_annotation_uniprot_helpers.R` |
| Does an undeclared extra change an identity? | No. `runDIAReferenceEnrichment` builds the identity and the product from `parameters[declared]` and the service only; the extras (`da_results`, `pathway_dir`, `go_annotations`, `process_fn`, and now an undeclared `go_propagation`) reach only the kernel through `runProjected` | `R/utils_workflow_dia_da.R`, `R/utils_workflow_operations.R` (`runProjected`) |
| Provenance class of `go_propagation` | `no_shipped_default` (no `config.ini` key; as `correction_method`); label "Gene Ontology sets include parent terms (relations followed)" (VAL-41) | `R/utils_workflow_operations.R` |
| Record path of `correction_method` | the enrichment step's declared parameter, `c("enrichment", "correction_method")`, read by `workflowRecordReportParameters` as a parameter; on the clusterProfiler route its value is the effective `"BH"` (A107-03), on the g:Profiler route the selection (`gSCS`, `fdr`) | `R/utils_workflow_analysis_record.R` |
| g:Profiler: correction per data source | __F2_PERSOURCE__ | |
| g:Profiler: `domain_scope = "custom"` | Confirmed. g:Profiler documentation (page `docs`, fetched in R0): "Custom over all known genes -- Use the set of all known genes in the custom background as the statistical domain", against "Custom over annotated genes -- Use the set of all genes that are annotated in the data source and also included in the custom background"; the gprofiler2 0.2.4 vignette maps `custom_bg` to `domain_scope = "custom"` and the annotated variant to `"custom_annotated"`. So the background is the custom background's identifiers that g:Profiler knows, annotated or not | `f1/gprofiler-page-docs.html.txt`; gprofiler2 vignette |
| The two citations | Confirmed (Europe PMC metadata, R0): Kolberg L, Raudvere U, Kuzmin I, Adler P, Vilo J, Peterson H. "g:Profiler-interoperable web service for functional enrichment analysis and gene identifier mapping (2023 update)". Nucleic Acids Research 2023, doi 10.1093/nar/gkad347 (the service). Kolberg L, Raudvere U, Kuzmin I, Vilo J, Peterson H. "gprofiler2 -- an R package for gene list functional enrichment analysis and namespace conversion toolset g:Profiler". F1000Research 2020, doi 10.12688/f1000research.24956.2 (the package; also the installed package's `CITATION`) | `f1/europepmc-*.json` |
| clusterProfiler 4.20.0 `enricher` on the lock | Confirmed as section 2.1 read it on devel: `enricher` builds a GSON object and calls `enrichit::ora_gson` (enrichit 0.2.0); defaults `minGSSize = 10`, `maxGSSize = 500`, `qvalueCutoff = 0.2`. In `ora_gson` the background is the universe intersected with the genes in TERM2GENE; set sizes are counted after that intersection; sets outside 10 to 500 are removed; rows with `Count > 0` kept; `p.adjust` (BH) and `qvalue` (`calculate_qvalue`) over those rows; `BgRatio` = set size / background size | `f1/enricher-on-lock.txt` |
