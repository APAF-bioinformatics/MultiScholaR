# KIT-167: a cut of the B36 GEO annotation

- `CUT-b36.gff`: 113 feature lines and the version line of `GSE152966_GCF_900622635.1_B36_genomic.gff.gz` (GEO series
  GSE152966; 365,810 B, SHA-256 `24d67f30592c2cbc6a0de9390eea3d856ff21ca089ac538aea125673aaee3024`; NCBI places no
  restrictions on the use or distribution of GEO data), the lines of 54 loci: the first 40 coding loci, CDS split over
  several rows, pseudogenes, RNA genes and protein accessions on several loci, as found.
- `CUT-b36-expected.tsv`: each locus's protein accessions, derived with awk from the cut's CDS rows (`locus_tag`,
  `protein_id`), not by the package's parser.
- Made on a runner by `docs/qa/workflow-execution/kit/kit-167/pass2-inputs.R` (section CUT), run 45386584323.
  SHA-256: `CUT-b36.gff` `7f9b5503fe1f903cdbc6a4fb883c74df4fa122563a1b1fcb07cbeb2f4f5914b3`, `CUT-b36-expected.tsv`
  `2462039a0f60b155e086a9a6f670c4984ce053cb210c6f53cf9ad17d1b9710dc`.
- The GEO GTF files of the four strains carry no `protein_id` (F4), so the parser refuses them by design; no GTF cut.
