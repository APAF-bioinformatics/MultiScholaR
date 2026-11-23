# KIT-173 fixtures: Gene Ontology parent terms on the proteomics UniProt route

Fetched once, on a GitHub-hosted runner, by `docs/qa/workflow-execution/kit/kit-173/f1-feasibility.R` (run R0b,
45121463951, commit 32ffcf04; the UniProt files are byte-identical to R0's, run 45119655641). Tests read them with the
network closed; `test-go-propagation.R` checks each file against the SHA-256 below.

| File | SHA-256 | Source and content |
|---|---|---|
| `uniprot-yeast-559292.tsv.gz` | `02a9ec5c8c51cab10a89c145c48f7addaff7395505c568338d16a9bf2d388d0e` | UniProt REST `uniprotkb/stream`, query `(organism_id:559292) AND (reviewed:true)`, fields `accession,go_id,gene_oln,reviewed`, TSV as served (SHA-256 d5e8ddda...067eaf), gzip; release 2026_04 (04-November-2026; `X-UniProt-Release`); 6,733 entries, 6,052 with GO ids |
| `uniprot-ecoli-83333.tsv.gz` | `e4800c1ca50ac9852e5ec8b40255f38336f4e601faa4b18b7771ee6f15d4cbe7` | as above for `(organism_id:83333) AND (reviewed:true)` (TSV SHA-256 704966f7...0a62e85); release 2026_04; 4,531 entries, 4,024 with GO ids |
| `gprofiler-scerevisiae-go.ENSG.gmt.gz` | `c1e06f30e2d60c94daaaaa8a298f8f885e58bd0b4fe49c2a2ac629b9ba17ca6b` | g:Profiler `https://biit.cs.ut.ee/gprofiler/static/gprofiler_full_scerevisiae.ENSG.gmt` (version e114_eg62_p19_9dc92ac8; GO classes releases/2026-01-23; file as served SHA-256 6153b378...460a37, 16,629 lines, sources GO, HP, WP): the 8,502 GO lines only, each line's description replaced by the term id, gzip. The default g:Profiler data, which include electronic GO annotations; no variant without them is served under the probed names (`f1/gprofiler-gmt-variants.tsv`). Genes are ordered locus names (6,589; none unmapped to the UniProt table) |
| `t7-annotation-ecoli.tsv` | `037e743cc19e831551210644ad24e6834eb1e8146e4b5eea27fbdd6094a83ae7` | 70 rows of the E. coli table (T7, T8, T3(b), T11(a)), selected by the F1 rule (`feasibility.md`): P = GO:0009067 (25 proteins after propagation, none direct), D = GO:0070180 (10 direct), 34 seeded fill proteins, and one root-only protein. **Edited:** P77326 has no GO annotation in UniProt 2026_04; here its `Gene.Ontology.IDs` is `GO:0008150` alone (no reviewed E. coli or yeast entry has only aspect-root annotations, so the root-only protein T7 and T8 need is synthetic) |
| `t7-fixture-entries.tsv` | `799bd4061718da7cf98e2c8a691f4ce5720530f15ddf6005ce60e0adae4a4a96` | the 70 entries and their role (`p_term`, `d_term`, `fill`, `root_only`) |
| `report-base/LFQ_report.R` | `56fe2b0e7b3983fbb5c44d9a52671052b9816abd57276e58919c2938ba053e5a` | the `GO Enrichment Section` chunk of `inst/reports/proteomics/LFQ_report.rmd` at a60c4e9a (base), verbatim, for T11(e) |
| `report-base/TMT_report.R` | `4c748776d7a67be5442079f46cc0c2945b4ccc92598355ec03f597472fd31864` | the same chunk of `TMT_report.rmd` at a60c4e9a |
| `report-base/DIANN_report.R` | `c2a74607f3d8580c461d223a659ef16e1d5601d1ba8d229c8348eb97bb765a28` | the same chunk of `DIANN_report.rmd` at a60c4e9a |
| `report-base/DIANN_limpa_report.R` | `4949e6a16db8650ca764c63d65ba8a4b48743361ad3f49e560bc319f43a73147` | the same chunk of `DIANN_limpa_report.rmd` at a60c4e9a |
| `replay-base/record.json` | `0812c62c1a6297bd0ed56524c1f4e60afc040ec2e2a6b5efe2a904700e08fff6` | T7(a) capture on the base code (below) |
| `replay-base/enricher.tsv` | `5df53ef7740f2a59bd9bb8e2cb5f6f15374caab8a1b9dee2c0fd0dd71bb442fd` | T7(a) capture on the base code (below) |
| `replay-base/evidence.rds` | `0e65813e03ecc12b34b0080408dccf0d7fc993008ce7fdc0f814b17afee16ce9` | T7(a) capture on the base code (below) |
| `replay-base/product.rds` | `65ab5efb635c3cb760612a658b4c9398bee51dc7a42727125a97734dc83ac7be` | T7(a) capture on the base code (below) |
| `replay-base/capture.tsv` | `5ebbb5539f9a18c17f89714b55ac43303a612bcdb535e5b906819aefafd31df3` | T7(a) capture on the base code (below) |

`replay-base/` (T7(a)): the outputs of `docs/qa/workflow-execution/kit/kit-173/capture-record.R` run on the
tests-only commit 536080c5 (R1, run of `k173-r1-536080c5`), whose `R/` equals the base a60c4e9a: `record.json`,
`enricher.tsv`, `evidence.rds`, `product.rds` and `capture.tsv` (commit, versions, declared names and the SHA-256 of each
file). The captured project store (`project.tar`) was withheld by the run's upload scan and is not committed; T7(b)
replays the base record on a project with the same import and design states built from the fixtures above.

Licences: UniProt CC BY 4.0 (The UniProt Consortium, UniProt: the Universal Protein Knowledgebase, Nucleic Acids
Res.); Gene Ontology annotations and g:Profiler's GO GMT, CC BY 4.0 (Gene Ontology Consortium; g:Profiler, Kolberg et
al. 2023, Nucleic Acids Res., doi 10.1093/nar/gkad347). The edit to P77326 above is ours.
