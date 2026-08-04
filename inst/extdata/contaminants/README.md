# Contaminant manifest policy

MultiScholaR does not currently redistribute a cRAP accession or sequence
manifest. The authoritative GPM cRAP page identifies the current list as
version `2012.01.01`, but it does not publish an explicit redistribution
licence for the downloadable FASTA/list.

Source: <https://www.thegpm.org/crap/>

Until source and redistribution terms are explicit, users can supply a
one-column, TSV, or data-frame accession manifest to
`classifyPeptideBiologicalExclusions()`. The audit output records the supplied
path/source and checksum. Explicit DIA-NN decoy/contaminant columns and
identifier tags remain supported without a manifest.

Do not add copied cRAP accessions or sequences to this directory without
recording source URL, version, retrieval date, checksum, and redistribution
licence.
