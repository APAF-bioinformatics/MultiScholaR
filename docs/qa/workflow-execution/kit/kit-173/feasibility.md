# KIT-173 feasibility (F1, F2)

Ticket revision 1.5. F1 ran on `run-tests` (R 4.6.1, reference BLAS, ubuntu-22.04) as R0 (run 45119655641, commit
36a34cab) and its correction pass R0b (commit 32ffcf04); the run ids, minutes and core-hours are in `runs.tsv`. The
script is `f1-feasibility.R`; its outputs are under `f1/` (R0b; R0's where only R0 has them, named `r0-`). F2 is a
reading of the code on the branch and of g:Profiler's published documentation as fetched by F1 (no call from the
laptop).

## F1: observed values

| Quantity | Observed (R0b; R0 identical where both ran) | Ticket (2.2, 2.3) |
|---|---|---|
| Versions on the lock | GO.db 3.23.1, AnnotationDbi 1.74.0, clusterProfiler 4.20.0, enrichit 0.2.0, gprofiler2 0.2.4 | as locked |
| `GO_dbInfo()` | go-basic.obo, https://current.geneontology.org/ontology/go-basic.obo, 2026-01-23 | 2026-01-23 |
| Relation labels the parent maps return (`toTable`, `as.list` names) | `isa`, `part of`, `regulates`, `positively regulates`, `negatively regulates` | the SQLite labels; the constants use these |
| Edges per label | BP isa 41,725, part of 4,493, regulates 2,867, positively 2,506, negatively 2,507; MF isa 12,413, part of 8; CC isa 4,664, part of 1,824 | same |
| Terms BP / MF / CC (+ `all`) | 24,547 / 10,123 / 4,069 (+1) | same |
| Term-ancestor pairs, five relations (recursive SQL) and `GO*ANCESTOR` | BP 464,934, MF 64,418, CC 42,075; `GO*ANCESTOR` identical pair by pair in each aspect | same |
| Pairs, `isa` + `part of` | BP 352,532, MF 64,418, CC 42,075 | same |
| Edges to another aspect | none; each root has one edge to `all` | same |
| GO:0006110 ancestors two / five; GO:0006096 and GO:0008152 among them | 18 (neither) / 75 (both); GO:0006096 has 57 under either | same |
| Secondary / obsolete ids; GO:0000002, GO:0000021 | 3,297 / 9,861; GO:0000002 obsolete; GO:0000021 secondary of GO:0007131 | same |
| `enricher` defaults | `minGSSize = 10`, `maxGSSize = 500` | same |
| UniProt fixtures | release 2026_04 (04-November-2026), yeast 6,733 entries (6,052 with GO), E. coli 4,531 (4,024 with GO); all reviewed; the same SHA-256 in R0 and R0b | fetched once |
| Root-only reviewed entries | yeast 0, E. coli 0 (no entry with GO:0008150 alone either) | "chosen in F1": none exists; see below |
| g:Profiler | version e114_eg62_p19_9dc92ac8 (Ensembl 114), GO classes releases/2026-01-23, annotations BioMart; GMT `gprofiler_full_scerevisiae.ENSG.gmt` (HTTP 200; `.name.gmt` 200; four probed no-IEA or zip names 404): 16,629 lines, 8,502 GO terms, all held by GO.db 3.23.1 | GMT with electronic annotations |
| Mapping by ordered locus name | 6,589 GO genes in the GMT, 0 unmapped (0 %); every yeast entry has an ordered locus name | loss reported; above 20 % a named cause |
| T5 pilot (50 terms of 2,703 in scope, seed 20261123; GO.db closures, not the function) | median C / P: five relations 0.937 / 0.989 (IQR C 0.808-1.000); two relations 0.879 / 1.000; none 0.017 / 1.000 (24 of 50 terms empty, P undefined) | rule: C, P >= 0.90; gain >= 0.20 |
| 9.1 planted BP terms (E. coli, shipped value, seed 20261123) | specific GO:0032790 ribosome disassembly (11 proteins, of 481 candidates); mid GO:0006396 RNA processing (127, of 175); broad GO:0044282 small molecule catabolic process (307, of 26) | one per band |
| T7 fixture (70 proteins; rule below) | P = GO:0009067 aspartate family amino acid biosynthetic process (0 direct, 25 propagated under either relation set); D = GO:0070180 large ribosomal subunit rRNA binding (10 direct); 87 terms reach 10 proteins only after propagation, 75 of them never annotated directly; root-only protein P77326 (synthetic, below) | >= 60 proteins; a parent term >= 10 only after propagation; a direct term >= 10; one root-only protein |

T7 selection rule (fixed in the script before any enrichment): P is drawn (seed 20261123 + 7) from the sorted BP
terms with at most 3 direct and 15 to 40 propagated proteins (five relations, roots removed) in the full E. coli
table; D is drawn next from the sorted terms with 10 to 20 direct proteins; the fixture is P's propagated proteins,
D's direct proteins, a seeded fill from the other proteins with a non-root term up to 69, and one root-only protein.

The pilot meets the rule's three conditions for the shipped value, so no named-cause revision of the rule is needed
(mapping loss 0 %, no reference term outside GO.db, and the GMT is the variant with electronic annotations).

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
| g:Profiler: correction per data source | Confirmed by the gprofiler2 0.2.4 vignette (the package's published documentation, as installed): the result metadata give "the statistics of data sources that are used in the enrichment testing", including for each source "the nominal significance threshold for this source" and its domain size, so each source is corrected on its own. The g:Profiler documentation page (R0) and the 2019 update paper's full text (R0b) describe the corrections without naming the family; the full texts of the 2023 and 2020 papers could not be read (Europe PMC answered HTTP 500 twice each in R0b; these two fetches were F2's documentation reads, not F1 fixtures, and are reported, not retried). The statement is printed on the strength of the vignette; the implementation audit may judge that insufficient, in which case the clause "each data source as its own family" is removed (F2's rule) | |
| g:Profiler: `domain_scope = "custom"` | Confirmed. g:Profiler documentation (page `docs`, fetched in R0): "Custom over all known genes -- Use the set of all known genes in the custom background as the statistical domain", against "Custom over annotated genes -- Use the set of all genes that are annotated in the data source and also included in the custom background"; the gprofiler2 0.2.4 vignette maps `custom_bg` to `domain_scope = "custom"` and the annotated variant to `"custom_annotated"`. So the background is the custom background's identifiers that g:Profiler knows, annotated or not | `f1/gprofiler-page-docs.html.txt`; gprofiler2 vignette |
| The two citations | Confirmed (Europe PMC metadata, R0): Kolberg L, Raudvere U, Kuzmin I, Adler P, Vilo J, Peterson H. "g:Profiler-interoperable web service for functional enrichment analysis and gene identifier mapping (2023 update)". Nucleic Acids Research 2023, doi 10.1093/nar/gkad347 (the service). Kolberg L, Raudvere U, Kuzmin I, Vilo J, Peterson H. "gprofiler2 -- an R package for gene list functional enrichment analysis and namespace conversion toolset g:Profiler". F1000Research 2020, doi 10.12688/f1000research.24956.2 (the package; also the installed package's `CITATION`) | `f1/europepmc-*.json` |
| clusterProfiler 4.20.0 `enricher` on the lock | Confirmed as section 2.1 read it on devel: `enricher` builds a GSON object and calls `enrichit::ora_gson` (enrichit 0.2.0); defaults `minGSSize = 10`, `maxGSSize = 500`, `qvalueCutoff = 0.2`. In `ora_gson` the background is the universe intersected with the genes in TERM2GENE; set sizes are counted after that intersection; sets outside 10 to 500 are removed; rows with `Count > 0` kept; `p.adjust` (BH) and `qvalue` (`calculate_qvalue`) over those rows; `BgRatio` = set size / background size | `f1/enricher-on-lock.txt` |
