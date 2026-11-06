# KIT-101: the offline UniProt cache for the entry-point suite's proteomics lanes. One row per accession of the e2e
# proteomics FASTA (tests/testdata/e2e/prot_dia/proteins.fasta: the protein universe every proteomics lane shares since
# revision 2.2, 400 synthetic accessions written by tests/testdata/e2e/generate_fixtures.R), in the columns
# `getUniprotAnnotations()` caches (`createEmptyUniprotTable()` plus the GO columns of the enrichment double). The design
# step's `getUniprotAnnotationsFull()` finds this file in `<results_dir>/cache/uniprot_annotations/` and does not query
# the network. Synthetic; deterministic (no random draws).
# Usage (repository root): Rscript tests/testdata/gui-entry/uniprot/build.R
headers <- grep("^>", readLines("tests/testdata/e2e/prot_dia/proteins.fasta"), value = TRUE)
entry <- sub("^>sp\\|([^|]+)\\|.*$", "\\1", headers)
gene <- sub("^.* GN=([^ ]+) .*$", "\\1", headers)
name <- sub("^>sp\\|[^|]+\\|[^ ]+ (.*) OS=.*$", "\\1", headers)
cache <- data.frame(Entry = entry, From = entry, gene_names = gene, Protein_existence = "Evidence at protein level",
  Protein_names = name, annotation_score = 5, Organism = "Homo sapiens",
  go_id_go_biological_process = "GO:0006915; GO:0007049", go_term_go_biological_process = "apoptotic process; cell cycle",
  go_id_go_molecular_function = "GO:0005515", go_term_go_molecular_function = "protein binding",
  go_id_go_cellular_compartment = "GO:0005739", go_term_go_cellular_compartment = "mitochondrion", stringsAsFactors = FALSE)
saveRDS(cache, "tests/testdata/gui-entry/uniprot/uniprot_annotations.RDS", version = 3)
cat("wrote", nrow(cache), "accessions\n")
