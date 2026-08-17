# Contaminant manifest policy

MultiScholaR does not currently redistribute a cRAP accession or sequence
manifest. The authoritative GPM cRAP page identifies the current list as
version `2012.01.01`, but it does not publish an explicit redistribution
licence for the downloadable FASTA/list.

Source: <https://www.thegpm.org/crap/>

Until source and redistribution terms are explicit, users can supply a TSV,
CSV, data-frame, or accession vector to
`classifyPeptideBiologicalExclusions()`. The current tabular manifest contract
uses these columns:

- `accession` (required)
- `manifest_schema_version` (currently `1.0.0`)
- `version` (the user resource version)
- `source_name` (a stable citation or resource name)
- `source_uri` and `license` (optional portable provenance)

The audit output records normalized accessions, contract/validation status,
source metadata, and deterministic SHA-256 fingerprints. It records only a
file basename as the input source and never persists an absolute local path.
Valid historical one-column files, data frames, and accession vectors remain
accepted through the explicit `legacy_adapter` contract. Missing path-like
inputs and malformed metadata fail explicitly. Explicit DIA-NN
decoy/contaminant columns and identifier tags remain supported without a
manifest.

The package test fixture contains only invented accessions and is not a copied
or derived third-party contaminant resource.

Do not add copied cRAP accessions or sequences to this directory without
recording source URL, version, retrieval date, checksum, and redistribution
licence.
